#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use Parallel::ForkManager;

# -----------------------------
# Usage / help message
# -----------------------------
sub usage {
    my ($msg) = @_;
    $msg ||= "";
    $msg .= "\n" if $msg ne "" && $msg !~ /\n$/;

    die $msg . <<"USAGE";
Usage:
  perl $0 --ais AIS_FILE --ref REFNAME --pairwisebam PAIR_FILE [--nproc N]

Required options:
  --ais          AIS file (SNP list) with header like:
                   Chr  POS  REF  ALT  Parent1  Parent2
                 and following columns per SNP.

  --ref          Name of a shell variable that holds the reference genome path.
                 For example:
                   export np7=/data/share/datesets/msu/nipponbare.fasta
                 then use:
                   --ref np7

  --pairwisebam  Tab-delimited text file with 3 columns:
                   <control_bam>   <pool_bam>   <output_file>
                 Example:
                   ctrl1.bam  pool1.bam  out1.txt
                   ctrl2.bam  pool2.bam  out2.txt

Optional:
  --nproc        Number of pair×chromosome jobs to run in parallel.
                 Default: 1 (no parallelism).

Example:
  export np7=/data/share/datesets/msu/nipponbare.fasta
  perl $0 \\
    --ais filtered_ais.txt \\
    --ref np7 \\
    --pairwisebam pairs.txt \\
    --nproc 4

USAGE
}

# -----------------------------
# Parse command line options
# -----------------------------
my ($ais_file, $refName, $pairwise_file);
my $nproc = 1;

GetOptions(
    'ais=s'         => \$ais_file,
    'ref=s'         => \$refName,
    'pairwisebam=s' => \$pairwise_file,
    'nproc=i'       => \$nproc,
) or usage("Error parsing command line options.");

# Check required options
my @missing;
push @missing, "--ais"         unless defined $ais_file;
push @missing, "--ref"         unless defined $refName;
push @missing, "--pairwisebam" unless defined $pairwise_file;

if (@missing) {
    usage("Missing required option(s): " . join(", ", @missing));
}

$nproc = 1 if !defined($nproc) || $nproc < 1;

# -----------------------------
# Load AIS into memory
#   - get parent names from header
#   - group sites by chromosome
# -----------------------------
my ($parent1_name, $parent2_name, $sites_by_chr_ref, $chr_order_ref) =
    load_ais($ais_file);

my %sites_by_chr = %{$sites_by_chr_ref};
my @chr_order    = @{$chr_order_ref};

# -----------------------------
# Global configuration
# -----------------------------
my $reference_genome = resolve_ref_path($refName);
my $samtools119      = "singularity exec -B /data/:/data/ /home/wangzicheng/work/0.script_software_docker/dockers/samtools119.sif samtools";

# -----------------------------
# Unique ID for temporary files
# -----------------------------
my $unique_id = time() . "_" . int(rand(10000));

# -----------------------------
# Build chr-specific region files
# -----------------------------
my %region_file_for_chr = create_region_files(\%sites_by_chr, \@chr_order, $unique_id);

# -----------------------------
# Read pairwise BAM list
#   Each line: <control_bam> <pool_bam> <output_file>
# -----------------------------
my @pairs;
open(my $pfh, '<', $pairwise_file) or die "Can't open $pairwise_file: $!\n";
while (<$pfh>) {
    chomp;
    next if /^\s*$/;
    next if /^\s*#/;

    my ($bam_control, $bam_pool, $out_file) = split;
    unless (defined $bam_control && defined $bam_pool && defined $out_file) {
        die "Invalid line $. in $pairwise_file: expect 3 columns: <control_bam> <pool_bam> <output_file>\n";
    }
    push @pairs, [$bam_control, $bam_pool, $out_file];
}
close $pfh;

# -----------------------------
# Parallel execution over (pair, chr)
# -----------------------------
my $pm = Parallel::ForkManager->new($nproc);

for (my $pair_idx = 0; $pair_idx <= $#pairs; $pair_idx++) {
    my ($bam_control, $bam_pool, $out_file) = @{$pairs[$pair_idx]};

    for my $chr (@chr_order) {
        my $region_file_chr = $region_file_for_chr{$chr};
        next unless defined $region_file_chr;

        $pm->start and next;  # child

        run_pair_chr_job(
            pair_idx      => $pair_idx,
            chr           => $chr,
            bam_control   => $bam_control,
            bam_pool      => $bam_pool,
            region_file   => $region_file_chr,
            sites_chr     => $sites_by_chr{$chr},
            uid           => $unique_id,
            samtools_cmd  => $samtools119,
            ref_genome    => $reference_genome,
        );

        $pm->finish;
    }
}

$pm->wait_all_children;

# -----------------------------
# Merge per-chr results into final outputs
# -----------------------------
for (my $pair_idx = 0; $pair_idx <= $#pairs; $pair_idx++) {
    my ($bam_control, $bam_pool, $out_file) = @{$pairs[$pair_idx]};

    open(my $outfh, '>', $out_file) or die "Can't write to $out_file: $!";

    print $outfh join(
        "\t",
        "Chr",
        "POS",
        "CTRL_" . $parent1_name,
        "CTRL_" . $parent2_name,
        "POOL_" . $parent1_name,
        "POOL_" . $parent2_name
    ), "\n";

    for my $chr (@chr_order) {
        my $tmp_res = tmp_result_file_name($unique_id, $pair_idx, $chr);
        next unless -e $tmp_res;

        open(my $rfh, '<', $tmp_res) or die "Can't open $tmp_res: $!";
        while (<$rfh>) {
            print $outfh $_;
        }
        close $rfh;

        unlink $tmp_res;
    }

    close $outfh;
    print "Finished: control=$bam_control, pool=$bam_pool -> $out_file\n";
}

# clean region files
for my $chr (@chr_order) {
    my $rf = $region_file_for_chr{$chr};
    unlink $rf if defined $rf && -e $rf;
}

print "All pairwise BAM comparisons completed.\n";

# ================== Subroutines ==================

sub resolve_ref_path {
    my ($refName) = @_;
    my $refPath = "\$" . $refName;
    $refPath = `echo $refPath`;
    chomp($refPath);

    if (!$refPath) {
        die "Failed to resolve reference genome path from ref name '$refName'. ".
            "Make sure a variable named '$refName' is defined in your shell, e.g.\n".
            "  export $refName=/path/to/reference.fasta\n";
    }
    return $refPath;
}

sub load_ais {
    my ($ais_file) = @_;

    open(my $fh, '<', $ais_file) or die "Can't open $ais_file: $!";
    my $header = <$fh>;
    die "AIS file $ais_file is empty\n" unless defined $header;
    chomp $header;
    my @h = split(/\s+/, $header);

    my $p1_name = $h[4] // "Parent1";
    my $p2_name = $h[5] // "Parent2";

    my %sites_by_chr;
    my @chr_order;
    my %seen_chr;

    while (<$fh>) {
        chomp;
        next if /^\s*$/;
        my ($chr, $pos, $ref, $alt, $p1_base, $p2_base) = split;
        next unless defined $chr && defined $pos && defined $ref && defined $p1_base && defined $p2_base;

        my $site = {
            chr           => $chr,
            pos           => $pos + 0,
            ref           => $ref,
            alt           => $alt,
            parent1_base  => $p1_base,
            parent2_base  => $p2_base,
        };

        push @{$sites_by_chr{$chr}}, $site;

        unless ($seen_chr{$chr}) {
            push @chr_order, $chr;
            $seen_chr{$chr} = 1;
        }
    }
    close $fh;

    for my $chr (keys %sites_by_chr) {
        my @sorted = sort { $a->{pos} <=> $b->{pos} } @{$sites_by_chr{$chr}};
        $sites_by_chr{$chr} = \@sorted;
    }

    return ($p1_name, $p2_name, \%sites_by_chr, \@chr_order);
}

sub create_region_files {
    my ($sites_by_chr_ref, $chr_order_ref, $uid) = @_;
    my %sites_by_chr = %{$sites_by_chr_ref};
    my @chr_order    = @{$chr_order_ref};

    my %region_file_for_chr;

    for my $chr (@chr_order) {
        my $sites = $sites_by_chr{$chr} or next;
        my $region_file = "region_${uid}_${chr}.txt";

        open(my $rfh, '>', $region_file) or die "Can't write to $region_file: $!";
        for my $s (@{$sites}) {
            my $pos = $s->{pos};
            print $rfh join("\t", $chr, $pos - 1, $pos), "\n";
        }
        close $rfh;

        $region_file_for_chr{$chr} = $region_file;
    }

    return %region_file_for_chr;
}

sub tmp_result_file_name {
    my ($uid, $pair_idx, $chr) = @_;
    return "tmp_${uid}_pair${pair_idx}_${chr}.txt";
}

sub run_pair_chr_job {
    my %args        = @_;
    my $pair_idx    = $args{pair_idx};
    my $chr         = $args{chr};
    my $bam_control = $args{bam_control};
    my $bam_pool    = $args{bam_pool};
    my $region_file = $args{region_file};
    my $sites_chr   = $args{sites_chr};
    my $uid         = $args{uid};
    my $samtools    = $args{samtools_cmd};
    my $ref_genome  = $args{ref_genome};

    my %mp_ctrl;
    my %mp_pool;

    my $cmd = "$samtools mpileup -f $ref_genome -l $region_file $bam_control $bam_pool";

    open(my $mp_fh, "-|", $cmd) or die "Failed to run mpileup for pair_idx=$pair_idx, chr=$chr: $!";

    while (<$mp_fh>) {
        chomp;
        next if /^\s*$/;
        my @f = split;
        next unless @f >= 9;

        my ($m_chr, $pos, $ref_base) = @f[0..2];
        my ($d1, $b1, $q1, $d2, $b2, $q2) = @f[3..8];

        my $key = "$m_chr:$pos";
        $mp_ctrl{$key} = { bases => $b1, ref => $ref_base, depth => $d1 };
        $mp_pool{$key} = { bases => $b2, ref => $ref_base, depth => $d2 };
    }
    close $mp_fh;

    my $results_chr = compute_allele_counts($sites_chr, \%mp_ctrl, \%mp_pool);

    my $tmp_res = tmp_result_file_name($uid, $pair_idx, $chr);
    open(my $out_fh, '>', $tmp_res) or die "Can't write to $tmp_res: $!";

    foreach my $r (@{$results_chr}) {
        printf $out_fh "%s\t%d\t%d\t%d\t%d\t%d\n",
            $r->{chr},
            $r->{pos},
            $r->{ctrl_p1},
            $r->{ctrl_p2},
            $r->{pool_p1},
            $r->{pool_p2};
    }

    close $out_fh;
}

sub count_alleles {
    my ($bases, $parent1_base, $parent2_base, $ref_base) = @_;
    my ($p1_count, $p2_count) = (0, 0);

    if ($parent1_base eq $ref_base) {
        $p1_count = () = $bases =~ /[.,]/g;
        $p2_count = () = $bases =~ /$parent2_base/gi;
    } else {
        $p1_count = () = $bases =~ /$parent1_base/gi;
        $p2_count = () = $bases =~ /[.,]/g;
    }

    return ($p1_count, $p2_count);
}

sub compute_allele_counts {
    my ($sites_chr, $mpileup_control, $mpileup_pool) = @_;
    my @results;

    for my $s (@{$sites_chr}) {
        my $chr          = $s->{chr};
        my $pos          = $s->{pos};
        my $ref          = $s->{ref};
        my $parent1_base = $s->{parent1_base};
        my $parent2_base = $s->{parent2_base};

        my $key = "$chr:$pos";

        my ($ctrl_p1, $ctrl_p2) = (0, 0);
        if (exists $mpileup_control->{$key}) {
            ($ctrl_p1, $ctrl_p2) = count_alleles(
                $mpileup_control->{$key}{bases},
                $parent1_base,
                $parent2_base,
                $ref
            );
        }

        my ($pool_p1, $pool_p2) = (0, 0);
        if (exists $mpileup_pool->{$key}) {
            ($pool_p1, $pool_p2) = count_alleles(
                $mpileup_pool->{$key}{bases},
                $parent1_base,
                $parent2_base,
                $ref
            );
        }

        push @results, {
            chr      => $chr,
            pos      => $pos,
            ctrl_p1  => $ctrl_p1,
            ctrl_p2  => $ctrl_p2,
            pool_p1  => $pool_p1,
            pool_p2  => $pool_p2,
        };
    }

    return \@results;
}

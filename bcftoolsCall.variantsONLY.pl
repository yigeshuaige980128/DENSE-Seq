#!/usr/bin/env perl
use strict;
use warnings;
use Parallel::ForkManager;
use File::Basename;
use Getopt::Long;

# -----------------------------
# Parse command line options
# -----------------------------
my ($bamList, $refName, $outPrefix, $nproc_opt);

GetOptions(
    "bamlist=s"   => \$bamList,
    "ref=s"       => \$refName,
    "outprefix=s" => \$outPrefix,
    "nproc=i"     => \$nproc_opt,
) or die "USAGE: perl $0 --bamlist <bam.list> --ref <refEnvName> --outprefix <prefix> --nproc <int>\n";

die "USAGE: perl $0 --bamlist <bam.list> --ref <refEnvName> --outprefix <prefix> --nproc <int>\n"
    unless defined $bamList && defined $refName && defined $outPrefix && defined $nproc_opt;

($nproc_opt =~ /^\d+$/ && $nproc_opt > 0)
    or die "[ERR] --nproc must be a positive integer\n";

# -----------------------------
# Tool paths
#   - bcftools via Singularity container
#   - bgzip/tabix from host PATH
# -----------------------------
my $sif      = "/data/wangzicheng/0.script_software_docker/dockers/bcftools.119.sif";
my $bcftools = "singularity exec -B /data/:/data/ $sif bcftools";
my $bgzip    = "bgzip";
my $tabix    = "tabix";

# Ensure bgzip/tabix are in host PATH
sub assert_in_path {
    my ($exe) = @_;
    my $p = `bash -lc 'command -v $exe 2>/dev/null'`;
    chomp($p);
    $p or die "[ERR] '$exe' not found in PATH. Please install it or add to PATH.\n";
    return $p;
}
my $bgzip_path = assert_in_path($bgzip);
my $tabix_path = assert_in_path($tabix);

# -----------------------------
# Resolve reference path from environment variable name
#   --ref gives the env var name, e.g. np7 -> $np7
# -----------------------------
my $refPathEnv = "\$" . $refName;
my $refPath    = `echo $refPathEnv`;
chomp($refPath);
(-s $refPath) or die "[ERR] Reference fasta not found: $refPath\n";

(-s "$refPath.fai")
  or die "[ERR] Missing FASTA index ($refPath.fai). Run: samtools faidx $refPath\n";

# -----------------------------
# Get chromosome list (in .fai order)
# -----------------------------
my @chroms;
open my $FAI, "$refPath.fai" or die $!;
while (<$FAI>) {
    next if /^\s*$/;
    my ($chr, $len) = (split)[0,1];
    next unless defined $chr && defined $len && $len > 0;
    push @chroms, $chr;
}
close $FAI;
@chroms or die "[ERR] No sequences found in $refPath.fai\n";

# -----------------------------
# Parallelism (number of processes for per-chromosome steps)
# -----------------------------
my $nproc = $nproc_opt;
$nproc = @chroms if $nproc > @chroms;   # Do not exceed chromosome count
warn "[INFO] Using $nproc parallel processes for per-chromosome jobs.\n";

# -----------------------------
# Output filename prefix
#   Final prefix: <outprefix>.<ref>
# -----------------------------
my $prefix = "$outPrefix.$refName";

# -----------------------------
# Per chromosome: mpileup -> call -> norm
#   Intermediate VCFs will be removed later after genome-wide concats
# -----------------------------
my $pm = Parallel::ForkManager->new($nproc);

CHR:
for my $chr (@chroms) {
    $pm->start and next CHR;

    my $mp   = "$prefix.$chr.mpileup.vcf.gz";
    my $call = "$prefix.$chr.mpileup.call.vcf.gz";
    my $norm = "$prefix.$chr.mpileup.call.norm.vcf.gz";

    # mpileup (add FORMAT fields such as AD/ADF/ADR)
    my $cmd1 = "$bcftools mpileup -d 1000 -Q 20 -q 30 " .
               "-a FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR " .
               "--bam-list $bamList --fasta-ref $refPath -r $chr -O z -o $mp";
    system($cmd1) == 0 or die "[CHR:$chr] mpileup failed\n";
    system("$tabix_path -f -p vcf $mp") == 0 or die "[CHR:$chr] mpileup index failed\n";

    # call
    my $cmd2 = "$bcftools call -m -v --annotate FORMAT/GQ -O z -o $call $mp";
    system($cmd2) == 0 or die "[CHR:$chr] call failed\n";
    system("$tabix_path -f -p vcf $call") == 0 or die "[CHR:$chr] call index failed\n";

    # norm
    my $cmd3 = "$bcftools norm -f $refPath -m +both -O z -o $norm $call";
    system($cmd3) == 0 or die "[CHR:$chr] norm failed\n";
    system("$tabix_path -f -p vcf $norm") == 0 or die "[CHR:$chr] norm index failed\n";

    $pm->finish(0);
}
$pm->wait_all_children;

# -----------------------------
# Collect per-chromosome VCFs for genome-wide concats
# -----------------------------
my (@mp_list, @call_list, @norm_list);
for my $chr (@chroms) {
    my $mp   = "$prefix.$chr.mpileup.vcf.gz";
    my $call = "$prefix.$chr.mpileup.call.vcf.gz";
    my $norm = "$prefix.$chr.mpileup.call.norm.vcf.gz";

    (-s $mp)   or die "[ERR] Missing per-chr mpileup VCF: $mp\n";
    (-s $call) or die "[ERR] Missing per-chr call VCF: $call\n";
    (-s $norm) or die "[ERR] Missing per-chr norm VCF: $norm\n";

    push @mp_list,   $mp;
    push @call_list, $call;
    push @norm_list, $norm;
}

# -----------------------------
# Genome-wide VCF filenames
# -----------------------------
my $all_mp   = "$prefix.allchr.mpileup.vcf.gz";
my $all_call = "$prefix.allchr.mpileup.call.vcf.gz";
my $all_norm = "$prefix.allchr.mpileup.call.norm.vcf.gz";

my $final_clean = "$prefix.allchr.clean.vcf.gz";
my $log         = "$prefix.allchr.vcf_filter.log";

# vcf_filter.pl path
my $vcf_filter = "/data/wangzicheng/0.script_software_docker/script/vcf_filter.pl";
(-s $vcf_filter) or die "[ERR] Missing vcf_filter.pl at $vcf_filter\n";

# -----------------------------
# Concat into genome-wide VCFs in parallel:
#   - mpileup        -> allchr.mpileup.vcf.gz
#   - call           -> allchr.mpileup.call.vcf.gz
#   - norm (+filter) -> allchr.mpileup.call.norm.vcf.gz -> allchr.clean.vcf.gz
# -----------------------------
my $pm_concat = Parallel::ForkManager->new(3);

for my $task (qw/mp call norm/) {
    $pm_concat->start and next;

    if ($task eq 'mp') {
        my $cmd_concat_mp = "$bcftools concat -a -O z -o $all_mp " . join(' ', @mp_list);
        system($cmd_concat_mp) == 0 or die "[ERR] bcftools concat (mpileup) failed\n";
        system("$tabix_path -f -p vcf $all_mp") == 0 or die "[ERR] allchr.mpileup index failed\n";
    }
    elsif ($task eq 'call') {
        my $cmd_concat_call = "$bcftools concat -a -O z -o $all_call " . join(' ', @call_list);
        system($cmd_concat_call) == 0 or die "[ERR] bcftools concat (call) failed\n";
        system("$tabix_path -f -p vcf $all_call") == 0 or die "[ERR] allchr.mpileup.call index failed\n";
    }
    elsif ($task eq 'norm') {
        my $cmd_concat_norm = "$bcftools concat -a -O z -o $all_norm " . join(' ', @norm_list);
        system($cmd_concat_norm) == 0 or die "[ERR] bcftools concat (norm) failed\n";
        system("$tabix_path -f -p vcf $all_norm") == 0 or die "[ERR] allchr.mpileup.call.norm index failed\n";

        # Once norm concat is done, immediately run vcf_filter
        my $cmd_final = qq{bash -lc "set -euo pipefail;
          perl $vcf_filter '$all_norm' '$bamList' 2>>'$log' | $bgzip_path -c > '$final_clean'
        "};
        system($cmd_final) == 0 or die "[ERR] vcf_filter.pl or compression failed, see $log\n";
    }

    $pm_concat->finish(0);
}
$pm_concat->wait_all_children;

# -----------------------------
# Remove per-chromosome intermediate files (VCF + index)
#   Keep only genome-wide allchr.* and final clean VCF
# -----------------------------
for my $chr (@chroms) {
    my $mp   = "$prefix.$chr.mpileup.vcf.gz";
    my $call = "$prefix.$chr.mpileup.call.vcf.gz";
    my $norm = "$prefix.$chr.mpileup.call.norm.vcf.gz";

    for my $f ($mp, "$mp.tbi", $call, "$call.tbi", $norm, "$norm.tbi") {
        unlink $f if defined $f && -e $f;
    }
}

print "[DONE] Final VCF: $final_clean\n";
print "[INFO] Kept genome-wide intermediates:\n";
print "       $all_mp\n";
print "       $all_call\n";
print "       $all_norm\n";

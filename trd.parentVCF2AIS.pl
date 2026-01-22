#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
use Getopt::Long;

# ================== Parse command-line options ==================
# Expected options:
#   --input   : VCF file (plain or gzipped)
#   --p1      : parent1 sample name (exactly as in VCF header)
#   --p2      : parent2 sample name (exactly as in VCF header)
my ($vcf_file, $parent1_name, $parent2_name);

GetOptions(
    'input=s' => \$vcf_file,
    'p1=s'    => \$parent1_name,
    'p2=s'    => \$parent2_name,
) or die "Usage: perl $0 --input <vcf_file> --p1 <parent1_sample_name> --p2 <parent2_sample_name>\n";

# Check required options
if ( !defined $vcf_file || !defined $parent1_name || !defined $parent2_name ) {
    die "Usage: perl $0 --input <vcf_file> --p1 <parent1_sample_name> --p2 <parent2_sample_name>\n";
}

# ================== Build output filename ==================
my $vcf_basename = fileparse($vcf_file, qr/\.[^\.]+$/);
my $output_file  = "${vcf_basename}.AIS.txt";

# ================== Open VCF (plain or gz) ==================
my $fh;

# If file ends with .gz, decompress on the fly with gzip -dc
if ($vcf_file =~ /\.gz$/) {
    open($fh, "-|", "gzip", "-dc", $vcf_file)
        or die "Can't gunzip $vcf_file: $!";
} else {
    open($fh, '<', $vcf_file)
        or die "Can't open $vcf_file: $!";
}

open(my $out_fh, '>', $output_file) or die "Can't write to $output_file: $!";

# ================== Find parent columns in #CHROM line ==================
my ($parent1_idx, $parent2_idx);
my $saw_header = 0;

while (my $line = <$fh>) {
    next if $line =~ /^##/;       # skip meta lines
    if ($line =~ /^#CHROM/) {
        $saw_header = 1;
        chomp $line;
        my @header_fields = split /\t/, $line;   # VCF is tab-separated

        # from column index 9 on are sample names
        for my $i (9 .. $#header_fields) {
            $parent1_idx = $i if $header_fields[$i] eq $parent1_name;
            $parent2_idx = $i if $header_fields[$i] eq $parent2_name;
        }

        last;
    }
}

die "ERROR: Cannot find #CHROM header line in VCF\n"
    unless $saw_header;

die "ERROR: Cannot find parent1 sample '$parent1_name' in VCF header\n"
    unless defined $parent1_idx;

die "ERROR: Cannot find parent2 sample '$parent2_name' in VCF header\n"
    unless defined $parent2_idx;

# ================== Print output header ==================
print $out_fh join(
    "\t",
    qw/Chr POS REF ALT/,
    $parent1_name,
    $parent2_name
), "\n";

# ================== Process variant lines ==================
while (<$fh>) {
    next if /^#/;   # skip any remaining header lines
    chomp;          # IMPORTANT: remove trailing newline so last sample column is clean

    my @fields = split /\t/;

    # First four columns are fixed in VCF:
    #   0: CHROM, 1: POS, 2: ID, 3: REF, 4: ALT, ...
    my ($chr, $pos, $ref, $alt) = @fields[0,1,3,4];

    # Raw genotype fields (may contain GT:AD:DP:...)
    my $parent1_field = $fields[$parent1_idx];
    my $parent2_field = $fields[$parent2_idx];

    # Extract GT from the genotype fields (take substring before the first ':')
    my ($parent1_gt) = split /:/, $parent1_field;
    my ($parent2_gt) = split /:/, $parent2_field;

    # Keep only SNPs where one parent is 0/0 and the other is 1/1
    if (($parent1_gt eq '0/0' && $parent2_gt eq '1/1') ||
        ($parent1_gt eq '1/1' && $parent2_gt eq '0/0')) {

        # Determine allele for parent1 and parent2
        my ($parent1_base, $parent2_base);
        if ($parent1_gt eq '0/0') {
            $parent1_base = $ref;
            $parent2_base = $alt;
        } else {
            $parent1_base = $alt;
            $parent2_base = $ref;
        }

        printf $out_fh "%s\t%d\t%s\t%s\t%s\t%s\n",
            $chr, $pos, $ref, $alt, $parent1_base, $parent2_base;
    }
}

close $fh;
close $out_fh;

print "VCF processing completed. Results saved to $output_file\n";

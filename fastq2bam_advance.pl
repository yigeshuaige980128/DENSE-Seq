#!/usr/bin/env perl
use strict;
use warnings;
use Parallel::ForkManager;
use Cwd;
use Getopt::Long;
use Pod::Usage;
use File::Spec;
use File::Basename qw(basename);

# -------------------------
# set options
# -------------------------
my $help  = 0;
my $clean = 0;
GetOptions('help|?' => \$help, 'clean' => \$clean) or pod2usage(2);
pod2usage(1) if $help;

# -------------------------
# parameter checking
# -------------------------
die "perl fastq2bam.pl <reads.list> <refname> <nthreads> <nfork>\n" unless @ARGV==4;
my ($readList, $refName, $cpu, $fork) = @ARGV;

# -------------------------
# set softwares
# -------------------------
my $fastqc0119     = "singularity exec -B /data/:/data/ /home/software/dockers/reseq.v2.sif fastqc";
my $trimmomatic039 = "singularity exec -B /data/:/data/ /home/wangzicheng/work/0.script_software_docker/dockers/trimmomatic039.sif trimmomatic";
my $bwa0717        = "singularity exec -B /data/:/data/ /home/software/dockers/reseq.v2.sif bwa";
my $samtools119    = "singularity exec -B /data/:/data/ /home/wangzicheng/work/0.script_software_docker/dockers/samtools119.sif samtools";
my $gatk381        = "singularity exec -B /data/:/data/ /home/wangzicheng/work/0.script_software_docker/dockers/gatk3.381.sif java -Djava.io.tmpdir=/data/wangzicheng/tmp -jar /usr/GenomeAnalysisTK.jar";
my $nchr = 12;
my $currentDir = getcwd;

# 用于避免 "@RG" 在双引号字符串中的插值
my $RG_TAG = '@RG';

# -------------------------
# forks
# -------------------------
my $fastqc1_fork     = Parallel::ForkManager->new($fork);
my $trimmomatic_fork = Parallel::ForkManager->new($fork);
my $fastqc2_fork     = Parallel::ForkManager->new($fork);
my $map_fork         = Parallel::ForkManager->new($fork);
# 合并阶段后面再 new 一个 $merge_fork

# -------------------------
# helper utilities & logs
# -------------------------
my $log_dir = "$currentDir/_runlogs";
mkdir $log_dir unless -d $log_dir;

sub _status_file    { return "$log_dir/$_[0].status"; }   # 只写“备注”一列内容
sub _fail_step_file { return "$log_dir/$_[0].fail";   }   # 记录首个失败步骤与原因（step \t msg）

# 只记录“第一个失败步骤”
sub _write_once_fail {
    my ($sample, $step, $msg) = @_;
    my $ff = _fail_step_file($sample);
    return if -e $ff;
    if (open my $FH, '>', $ff) {
        print $FH "$step\t$msg\n";
        close $FH;
    } else {
        warn "[WARN] Cannot write $ff: $!";
    }
}

# 样品最终状态（第二列表达式），例如：
#   OK: /abs/path/to/sample.ref.highQual.bam
#   ERROR: merge_bams failed (...)
sub _note_status {
    my ($sample, $remark) = @_;
    my $sf = _status_file($sample);
    if (open my $FH, '>', $sf) {
        print $FH "$remark\n";
        close $FH;
    } else {
        warn "[WARN] Cannot write $sf: $!";
    }
}

sub _abs_path {
    my ($path) = @_;
    return File::Spec->file_name_is_absolute($path) ? $path : File::Spec->catfile($currentDir, $path);
}

# 运行命令，失败时记录并返回 0；成功返回 1
sub run_step {
    my (%opt) = @_;
    my $sample = $opt{sample} // 'NA';
    my $step   = $opt{step}   // 'unknown_step';
    my $cmd    = $opt{cmd}    // '';
    my $onfail = $opt{onfail} // '';

    if (!$cmd) {
        _write_once_fail($sample, $step, "empty command");
        return 0;
    }
    my $ret = system($cmd);
    if ($ret != 0) {
        my $msg = "$step failed (exit=$ret). $onfail";
        _write_once_fail($sample, $step, $msg);
        warn "[ERROR][$sample] $msg\n";
        return 0;
    }
    return 1;
}

# -------------------------
# parse <refname> 变量
# -------------------------
my $refPath = "\$".$refName;
$refPath = `echo $refPath`;
chomp ($refPath);

# -------------------------
# parse <reads.list>
#   %RG: fastqName -> sampleName
#   %RGreverse: sampleName -> [fastqName...]
# -------------------------
my @fastqNameFullPath;
my @fastqName;
my %RG;
my %RGreverse;
open (my $fh, $readList) or die "cannot open $readList\n";
while (<$fh>){
    chomp;
    my ($fastqNameFullPath, $sampleName) = split;
    next unless defined $fastqNameFullPath && defined $sampleName;
    push @fastqNameFullPath, $fastqNameFullPath;
    my ($fastqName) = $fastqNameFullPath =~ m|([^/]+)$|;
    push @fastqName, $fastqName;
    $RG{$fastqName} = $sampleName;
}
close $fh;

while (my ($fq, $sampleName) = each %RG) {
    $RGreverse{$sampleName} = [] unless exists $RGreverse{$sampleName};
    push @{$RGreverse{$sampleName}}, $fq;
}

foreach my $sampleName (keys %RGreverse) {
    print "sample: $sampleName, fastq files: @{$RGreverse{$sampleName}}\n";
}

# -------------------------
# create symbolic link
# -------------------------
for (@fastqNameFullPath) {
    my $fastq_mate1_pattern = "$_\_*1*";
    my $fastq_mate2_pattern = "$_\_*2*";
    my @fastq_mate1_file = glob($fastq_mate1_pattern);
    my @fastq_mate2_file = glob($fastq_mate2_pattern);
    my ($fastqName) = $_ =~ m|([^/]+)$|;
    print "creating symbolic links for $fastqName:\n";

    if (@fastq_mate1_file == 0) {
        warn "No files found for pattern $fastq_mate1_pattern\n";
    } else {
        print "Found mate1 files: @fastq_mate1_file\n";
    }

    if (@fastq_mate2_file == 0) {
        warn "No files found for pattern $fastq_mate2_pattern\n";
    } else {
        print "Found mate2 files: @fastq_mate2_file\n";
    }

    for my $fastq_mate1 (@fastq_mate1_file) {
        my $result = system("ln -s $fastq_mate1");
        warn "Failed to create symlink for $fastq_mate1: $!" if $result != 0;
    }

    for my $fastq_mate2 (@fastq_mate2_file) {
        my $result = system("ln -s $fastq_mate2");
        warn "Failed to create symlink for $fastq_mate2: $!" if $result != 0;
    }
}

# -------------------------
# trim raw fastq  (并行，失败不影响其他样品)
# -------------------------
for (my $i=0; $i<@fastqName; $i++) {
    $trimmomatic_fork->start and next;
    my $currentFastqName = $fastqName[$i];
    my $sample = $RG{$currentFastqName} // $currentFastqName;

    my $fastq_mate1_pattern = "${currentFastqName}_*1*";
    my $fastq_mate2_pattern = "${currentFastqName}_*2*";
    my $fastq_mate1_file = (glob($fastq_mate1_pattern))[0] // '';
    my $fastq_mate2_file = (glob($fastq_mate2_pattern))[0] // '';

    if (!$fastq_mate1_file || !$fastq_mate2_file) {
        _write_once_fail($sample, "trimmomatic", "missing input FASTQ for $currentFastqName");
        $trimmomatic_fork->finish; next;
    }

    my $cmd = "$trimmomatic039 PE -threads $cpu -phred33 $fastq_mate1_file $fastq_mate2_file ${currentFastqName}_1.paired.clean.fq.gz ${currentFastqName}_1.unpaired.clean.fq.gz ${currentFastqName}_2.paired.clean.fq.gz ${currentFastqName}_2.unpaired.clean.fq.gz ILLUMINACLIP:/home/wangzicheng/work/0.script_software_docker/trimmomatic039_adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36";

    run_step(sample=>$sample, step=>"trimmomatic", cmd=>$cmd)
        or do { $trimmomatic_fork->finish; next; };

    unlink("${currentFastqName}_1.unpaired.clean.fq.gz");
    unlink("${currentFastqName}_2.unpaired.clean.fq.gz");
    $trimmomatic_fork->finish;
}
$trimmomatic_fork->wait_all_children;

# -------------------------
# quality check (不等待——按你的要求)
# -------------------------
for (my $i=0; $i<@fastqName; $i++) {
    $fastqc2_fork->start and next;
    my $currentFastqName = $fastqName[$i];
    my $fastq_mate1_file = "${currentFastqName}_1.paired.clean.fq.gz";
    my $fastq_mate2_file = "${currentFastqName}_2.paired.clean.fq.gz";
    my $outputDir = $currentDir."/fastqCleanQual";
    mkdir $outputDir unless -d $outputDir;
    system "$fastqc0119 -t $cpu -o $outputDir $fastq_mate1_file $fastq_mate2_file";
    $fastqc2_fork->finish;
}
# 不 wait_all_children —— 按你的意思

# -------------------------
# map each fastq file; @RG in bam is fastq file name  (并行 + 稳健)
# -------------------------
for (my $i=0; $i<@fastqName; $i++) {
    $map_fork->start and next;
    my $currentFastqName = $fastqName[$i];
    my $sample = $RG{$currentFastqName} // $currentFastqName;

    my $fastq_mate1_file = "${currentFastqName}_1.paired.clean.fq.gz";
    my $fastq_mate2_file = "${currentFastqName}_2.paired.clean.fq.gz";
    my $RGsampleName = $currentFastqName.".$refName";

    # 检查输入
    unless (-e $fastq_mate1_file && -e $fastq_mate2_file) {
        _write_once_fail($sample, "bwa_mem", "trimmed FASTQ missing for $currentFastqName");
        $map_fork->finish; next;
    }

    # 1) bwa mapping  —— 这里用 $RG_TAG 拼接，避免 "@RG" 插值
    my $bwa_cmd = $bwa0717
                . " mem -t $cpu -M -v 1 -R '"
                . $RG_TAG
                . "\\tID:$RGsampleName\\tLB:$RGsampleName\\tPL:MGI\\tPU:$RGsampleName\\tSM:$RGsampleName' "
                . "$refPath $fastq_mate1_file $fastq_mate2_file | $samtools119 view -@ $cpu -Sb | "
                . "$samtools119 sort -@ $cpu -o $currentFastqName.$refName.bam -";
    my $index_bam_cmd = "$samtools119 index $currentFastqName.$refName.bam";
    run_step(sample=>$sample, step=>"bwa_mem", cmd=>"$bwa_cmd && $index_bam_cmd",
             onfail=>"FASTQ=$currentFastqName") or do { $map_fork->finish; next; };

    unlink($fastq_mate1_file);
    unlink($fastq_mate2_file);

    # 2) samtools properPair
    my $properPair_cmd = "$samtools119 view -@ $cpu -h -f 2 -b $currentFastqName.$refName.bam | $samtools119 sort -@ $cpu -o $currentFastqName.$refName.properPair.bam -";
    my $index_properPair_bam_cmd = "$samtools119 index $currentFastqName.$refName.properPair.bam";
    run_step(sample=>$sample, step=>"proper_pair", cmd=>"$properPair_cmd && $index_properPair_bam_cmd")
        or do { $map_fork->finish; next; };

    if ($clean) {
        unlink("$currentFastqName.$refName.bam");
        unlink("$currentFastqName.$refName.bam.bai");
    }

    # 3) multiple quality control
    my $multiqc_cmd = "$samtools119 view -@ $cpu -h -q 30 -F 4 -b $currentFastqName.$refName.properPair.bam | $samtools119 view -@ $cpu -h -F 256 -b - | $samtools119 view -@ $cpu -h -F 512 -b - | $samtools119 view -@ $cpu -h -F 1024 -b - | $samtools119 view -@ $cpu -h -F 2048 -b - | $samtools119 sort -@ $cpu -o $currentFastqName.$refName.properPair.multiqc.bam -";
    my $index_properPair_multiqc_bam_cmd = "$samtools119 index $currentFastqName.$refName.properPair.multiqc.bam";
    run_step(sample=>$sample, step=>"multi_qc", cmd=>"$multiqc_cmd && $index_properPair_multiqc_bam_cmd")
        or do { $map_fork->finish; next; };

    if ($clean) {
        unlink("$currentFastqName.$refName.properPair.bam");
        unlink("$currentFastqName.$refName.properPair.bam.bai");
    }

    # 4) remove soft clip (外部 perl 脚本)
    my $sam_get_soft_clip_readname_cmd = "$samtools119 view -h $currentFastqName.$refName.properPair.multiqc.bam | perl /home/wangzicheng/work/0.script_software_docker/script/sam_get_soft_clip_readname.pl - > $currentFastqName.$refName.softclip.rname";
    my $sam_mask_with_read_name = "$samtools119 view -h $currentFastqName.$refName.properPair.multiqc.bam | perl /home/wangzicheng/work/0.script_software_docker/script/sam_mask_with_read_name.pl - $currentFastqName.$refName.softclip.rname | $samtools119 view -@ $cpu -Sb - | $samtools119 sort -@ $cpu -o $currentFastqName.$refName.properPair.multiqc.rmSC.bam -";
    my $index_properPair_multiqc_rmSC_bam_cmd = "$samtools119 index $currentFastqName.$refName.properPair.multiqc.rmSC.bam";
    run_step(sample=>$sample, step=>"remove_softclip",
             cmd=>"$sam_get_soft_clip_readname_cmd && $sam_mask_with_read_name && $index_properPair_multiqc_rmSC_bam_cmd")
        or do { $map_fork->finish; next; };

    unlink("$currentFastqName.$refName.softclip.rname");

    if ($clean) {
        unlink("$currentFastqName.$refName.properPair.multiqc.bam");
        unlink("$currentFastqName.$refName.properPair.multiqc.bam.bai");
    }

    # 5) mark & remove duplicates
    my $sortName_cmd = "$samtools119 sort -n -@ $cpu -O BAM $currentFastqName.$refName.properPair.multiqc.rmSC.bam";
    my $fixmate_cmd  = "$samtools119 fixmate -m - -";
    my $sortPos_cmd  = "$samtools119 sort -@ $cpu -O BAM -";
    my $markdup_cmd  = "$samtools119 markdup -r -@ $cpu - - | $samtools119 sort -@ $cpu -o $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.bam -";
    my $index_properPair_multiqc_rmSC_rmDup_bam_cmd = "$samtools119 index $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.bam";
    run_step(sample=>$sample, step=>"markdup",
             cmd=>"$sortName_cmd | $fixmate_cmd | $sortPos_cmd | $markdup_cmd && $index_properPair_multiqc_rmSC_rmDup_bam_cmd")
        or do { $map_fork->finish; next; };

    if ($clean) {
        unlink("$currentFastqName.$refName.properPair.multiqc.rmSC.bam");
        unlink("$currentFastqName.$refName.properPair.multiqc.rmSC.bam.bai");
    }

    # 6) GATK3 indel realignment
    my $indel_realignment_cmd = "$gatk381 -T RealignerTargetCreator -R $refPath -I $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.bam -o $currentFastqName.$refName.IndelRealigner.intervals && $gatk381 -T IndelRealigner -R $refPath -I $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.bam -o $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.realign.tmp.bam --targetIntervals $currentFastqName.$refName.IndelRealigner.intervals && $samtools119 sort -@ $cpu -o $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.realign.bam $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.realign.tmp.bam";
    my $index_properPair_multiqc_rmSC_rmDup_realign_bam_cmd = "$samtools119 index $currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.realign.bam";
    run_step(sample=>$sample, step=>"indel_realign",
             cmd=>"$indel_realignment_cmd && $index_properPair_multiqc_rmSC_rmDup_realign_bam_cmd")
        or do { $map_fork->finish; next; };

    unlink("$currentFastqName.$refName.IndelRealigner.intervals");
    unlink("$currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.realign.tmp.bai");
    unlink("$currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.realign.tmp.bam");

    if ($clean) {
        unlink("$currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.bam");
        unlink("$currentFastqName.$refName.properPair.multiqc.rmSC.rmDup.bam.bai");
    }

    $map_fork->finish;
}
$map_fork->wait_all_children;

# -------------------------
# merge those bam & modify @RG (并行)
# rename those samples do not need merge
# -------------------------
my $merge_fork = Parallel::ForkManager->new($fork);

for my $sampleName (keys %RGreverse) {
    $merge_fork->start and next;

    my @expected = map { "$_.$refName.properPair.multiqc.rmSC.rmDup.realign.bam" } @{$RGreverse{$sampleName}};
    my @exists   = grep { -e $_ } @expected;

    if (@exists == 0) {
        # 前序流程没有任何 realign BAM
        my $note = "no realign BAMs generated";
        # 尝试附带首个失败步骤信息
        my $ff = _fail_step_file($sampleName);
        if (-e $ff) {
            if (open my $F, '<', $ff) {
                my $line = <$F>;
                close $F;
                chomp $line if defined $line;
                $note .= "; first failure: $line" if defined $line;
            }
        }
        _note_status($sampleName, "ERROR: $note");
        $merge_fork->finish; next;
    }

    my $final_bam = "${sampleName}.$refName.highQual.bam";

    if (@exists > 1) {
        # 需要合并（仅合并真正存在的）
        my $merged = "${sampleName}_merge.$refName.properPair.multiqc.rmSC.rmDup.realign.bam";
        my $list   = join(' ', @exists);

        run_step(sample=>$sampleName, step=>"merge_bams",
                 cmd=>"$samtools119 merge $merged $list",
                 onfail=>"Missing inputs: ".join(',', grep { !-e $_ } @expected))
            or do { _note_status($sampleName, "ERROR: merge_bams failed"); $merge_fork->finish; next; };

        # 这里同样用 $RG_TAG，避免 "@RG" 插值
        run_step(sample=>$sampleName, step=>"addreplacerg",
                 cmd=>"$samtools119 addreplacerg -m overwrite_all -r '"
                      . $RG_TAG
                      . "\\tID:$sampleName.$refName\\tLB:$sampleName.$refName\\tPL:MGI\\tPU:$sampleName.$refName\\tSM:$sampleName.$refName' -o $final_bam $merged")
            or do { _note_status($sampleName, "ERROR: addreplacerg failed"); $merge_fork->finish; next; };

        run_step(sample=>$sampleName, step=>"index_final",
                 cmd=>"$samtools119 index $final_bam")
            or do { _note_status($sampleName, "ERROR: index_final failed"); $merge_fork->finish; next; };

        unlink $merged;

        if ($clean) {
            for my $f (@exists) {
                unlink $f;
                unlink("$f.bai");
            }
        }
    }
    else {
        # 单个 BAM：复制 + addreplacerg
        my $one   = $exists[0];
        my $cpTMP = "${sampleName}_cp.$refName.properPair.multiqc.rmSC.rmDup.realign.bam";
        system "cp $one $cpTMP";  # 复制失败也让后一步去报错

        run_step(sample=>$sampleName, step=>"addreplacerg",
                 cmd=>"$samtools119 addreplacerg -m overwrite_all -w -r '"
                      . $RG_TAG
                      . "\\tID:$sampleName.$refName\\tLB:$sampleName.$refName\\tPL:MGI\\tPU:$sampleName.$refName\\tSM:$sampleName.$refName' -o $final_bam $cpTMP")
            or do { _note_status($sampleName, "ERROR: addreplacerg failed"); $merge_fork->finish; next; };

        run_step(sample=>$sampleName, step=>"index_final",
                 cmd=>"$samtools119 index $final_bam")
            or do { _note_status($sampleName, "ERROR: index_final failed"); $merge_fork->finish; next; };

        unlink $cpTMP;

        if ($clean) {
            unlink $one;
            unlink "$one.bai";
        }
    }

    # 成功：记录绝对路径
    my $abs_final = _abs_path($final_bam);
    _note_status($sampleName, "OK: $abs_final");

    $merge_fork->finish;
}
$merge_fork->wait_all_children;

# ----------------------------
# 汇总写出 fastq2bam.summary.tsv  （两列：Sample \t Remark）
# ----------------------------
my $summary = File::Spec->catfile($currentDir, "fastq2bam.summary.tsv");
if (open my $SUM, '>', $summary) {
    print $SUM "Sample\tRemark\n";
    for my $sampleName (sort keys %RGreverse) {
        my $remark;

        # 优先读取 status 文件（如果合并阶段明确写了）
        my $sf = _status_file($sampleName);
        if (-e $sf) {
            if (open my $F, '<', $sf) {
                my $line = <$F>;
                close $F;
                chomp $line if defined $line;
                $remark = $line if defined $line;
            }
        }

        # 没有 status，则根据最终 bam 是否存在推断
        if (!defined $remark) {
            my $final_bam = _abs_path("${sampleName}.$refName.highQual.bam");
            if (-e $final_bam) {
                $remark = "OK: $final_bam";
            } else {
                # 再退一步读取 fail_step
                my $ff = _fail_step_file($sampleName);
                if (-e $ff) {
                    if (open my $F, '<', $ff) {
                        my $line = <$F>;
                        close $F;
                        chomp $line if defined $line;
                        if (defined $line && $line =~ /^([^\t]+)\t(.*)$/) {
                            $remark = "ERROR at $1: $2";
                        }
                    }
                }
            }
        }

        $remark = "ERROR: unknown (no outputs)" unless defined $remark;
        print $SUM "$sampleName\t$remark\n";
    }
    close $SUM;
    print STDERR "[INFO] Summary written: $summary\n";
} else {
    warn "[WARN] Cannot write $summary: $!";
}

__END__

=pod

=head1 NAME

fastq2bam.pl

=head1 SYNOPSIS

 fastq2bam.pl [options] <reads.list> <refname> <nthreads> <nfork>

 Note: There should be no extra blank lines in the read.list file
 Note: The set number of <nfork> should not exceed the number of fastq file

 **the script will merge FASTQ files belonging to same sample
 **output includes BAM files for each sample and A folder displaying the quality of trimmed fastq data

 **fastq file trim criteria:
      1.sliding window trimming, SLIDINGWINDOW: 4:15
      2.trimming Low-Quality reads at the front and back ends, LEADING: 3, TRAILING: 3
      3.adapter contamination, adapter - 2 is seedMismatches, 30 is palindromeClipThreshold, the palindrome clip threshold. 10 is simpleClipThreshold, the simple clip threshold
      4.minimum length filtering, MINLEN: 36

 **bam file filtering criteria:
      1.mapping quality score below 30
      2.reads not properly paired (FLAG:2)
      3.reads unmapped (FLAG:4)
      4.reads marked as secondary alignment (FLAG:256)
      5.reads map to multiple locations (FLAG:2048)
      6.reads not passing filters (FLAG:512, placeholder)
      7.reads marked as duplicates (FLAG:1024, placeholder)
      8.reads has total soft-clipp segment greater than 10 (perl script)

 Options:
    -help               print the manual
    -clean              remove all intermediate files

 <reads.list> format (no header), first column list fastq name (accept full path), second column list sample name:

    sample1                         sample1
    sample2                         sample2
    sample3                         sample3
    sample4_1                       sample4 (need merge)
    sample4_2                       sample4 (need merge)
    sample4_3                       sample4 (need merge)
    sample5_1                       sample5 (need merge)
    sample5_2                       sample5 (need merge)

 <refname> requires a variable name that points to an indexed genome fasta file:
    eg. "np7=/data/share/datesets/msu/nipponbare.fasta"

 <nthreads> set how many thread can use, trimmomatic, bwa and samtools need

 <nfork> set the number of parallel task

=head1 CHANGES

- Fix: avoid unintended interpolation of @RG by composing strings with $RG_TAG.
- Add robust error handling: any step failure only marks that sample and continues.
- Parallelize merging / @RG rewriting per sample.
- Produce two-column log: fastq2bam.summary.tsv (Sample, Remark).

=head1 AUTHOR

wangGangDan (original)
Edits for robustness & parallel merge by ChatGPT

=cut

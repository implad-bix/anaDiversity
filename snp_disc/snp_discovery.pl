#!/usr/bin/env perl

use FindBin qw( $Bin );
use strict;

my $project_id  = shift;
my $working_dir = shift;
my $fasta_in    = shift;
my $fastq_in    = shift;
my $cmd = "";

print STDERR "\nExample run: perl snp_discovery.pl 1234 /tmp/o30 sample_files/sample.fas sample_files/sample_fq1.tar.gz\n";

die "\n\nUsage: perl $0 working_dir project_id file_ref.fasta file_fastq\n\n" unless ($project_id && -e $fasta_in && (-e $fastq_in || $fastq_in =~ /,/));
print STDERR "$project_id\t$fasta_in\t$fastq_in\n";
###1. The reads were mapped to the reference using BWA.
###(1) Index the reference with BWA and samtools, respectively.

#my $working_dir = "/tmp/dir_".$project_id;
my $samtools = $Bin."/exe/samtools-0.1.18/samtools";
my $bwa      = $Bin."/exe/bwa-0.7.17/bwa";
my $bcftools = $Bin."/exe/samtools-0.1.18/bcftools/bcftools"; 

$working_dir =~ s/\/+$//;
`mkdir $working_dir` unless (-e $working_dir);
my $ref      = $working_dir."/".$project_id.".fas";
my $fastq    = $working_dir."/".$project_id.".fq";
`cp $fasta_in $ref` unless (-e $ref);
#`cp $fastq_in $fastq` unless (-e $fastq);

`$bwa index -a is $ref`; 
#`bwa index -a is $ref`; 
`$samtools faidx $ref`;

#bwa index -a is reference.fasta
#samtools faidx reference.fasta
###(2)The reads were mapped to the reference using BWA.

print STDERR "\n\nStart bwa mem\n\n";

my $fastq = $fastq_in;
$fastq =~ s/,/ /g;
$cmd = "$bwa mem $ref $fastq > $working_dir/aln.sam";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.sam");

###(3)将sam格式的文件转换成bam格式的文件.
print STDERR "\n\nStart samtool view\n\n";

$cmd = "$samtools view -bS $working_dir/aln.sam > $working_dir/aln.bam";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.bam");
 
# -b output BAM
# -S input SAM

###(4)对bam文件进行排序
print STDERR "\n\nStart samtools sort\n\n";
$cmd = "$samtools sort $working_dir/aln.bam $working_dir/aln.sorted";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.sorted");

###(5)对排序后的bam文件，建立索引，生成后缀为.bai的文件，用于快速的随机处理。
print STDERR  "\n\nStart samtools index\n\n";
$cmd = "$samtools index $working_dir/aln.sorted.bam";
print STDERR $cmd, "\n\n";
`$cmd`;

###(6)Remove potential PCR duplicates	
print STDERR "\n\nStart samtools rmdup\n\n";
$cmd = "$samtools rmdup $working_dir/aln.sorted.bam $working_dir/aln.sorted.rmdup.bam";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.sorted.rmdup.bam");

`$samtools index $working_dir/aln.sorted.rmdup.bam`;

###2. 将bam文件生成bcf格式文件。
print STDERR "\n\nStart samtools mpileup\n\n";
$cmd = "$samtools mpileup -guSDf $ref $working_dir/aln.sorted.rmdup.bam > $working_dir/aln.bcf";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.bcf");

###3. 调用bcftools进行SNP和INDELpcalling
print STDERR "\n\nStart bcftools\n\n";
$cmd = "$bcftools view $working_dir/aln.bcf > $working_dir/snp_indel.vcf";
#$cmd = $Bin."/../samtools-0.1.18/bcftools/bcftools view $working_dir/aln.bcf > $working_dir/snp_indel.vcf";
#$cmd = "/home/cliu/my_apps/samtools-0.1.18/bcftools/bcftools view  $working_dir/aln.bcf > $working_dir/snp_indel.vcf";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/snp_indel.vcf");

###4. 对variant calling 的结果进行过滤
print STDERR "\n\nCompleted!\n\n";
#print STDERR "\n\nFilter results\n\n";
#`perl -ne 'print $_ if /DP4=(\d+),(\d+),(\d+),(\d+)/ && ($3+$4)>=10 && ($3+$4)/($1+$2+$3+$4)>=0.1' $working_dir/snp_indel.vcf > $working_dir/snp_indel.final.vcf`;

# DP4提供了四个数据：1 和正链一致的reads数，2 和负链一致的reads数，3 和正链上的variant上的reads数，4 和负链上的variant一致的reads数。
###5. 根据snp以及indel位点设计引物，验证SNP
#`nvcf2primer3in.pl $working_dir/snp_indel.final.vcf $ref > $working_dir/snp_indel.p3in`;
#`primer3_core -output=$working_dir./snp_indel.p3out $working_dir./snp_indel.p3in`;


#!/usr/bin/env perl

use FindBin qw( $Bin );
use strict;

my $project_id  = shift;
my $working_dir = shift;
my $fasta_in    = shift;
my $fastq_in    = shift;
my $cmd = "";


print STDERR "
Copyright (C) 2019 Linchun Shi and Chang Liu
           Institute of Medicinal Plant Development
	   Chinese Academy of Medical Science
Freely distributed under the GNU General Public License (GPLv3)
";

print STDERR "\nExample run: perl snp_discovery.pl 1234 /tmp/o30 sample_files/sample.fas sample_files/sample_fq1.tar.gz\n";

die "\nUsage: perl $0 working_dir project_id file_ref.fasta file_fastq\n\n" unless ($project_id && -e $fasta_in && (-e $fastq_in || $fastq_in =~ /,/));
print STDERR "$project_id\t$fasta_in\t$fastq_in\n";

###(1) set up the working environment
my $samtools = $Bin."/exe/samtools-0.1.18/samtools";
my $bwa      = $Bin."/exe/bwa-0.7.17/bwa";
my $bcftools = $Bin."/exe/samtools-0.1.18/bcftools/bcftools"; 

$working_dir =~ s/\/+$//;
`mkdir $working_dir` unless (-e $working_dir);
my $ref      = $working_dir."/".$project_id.".fas";
my $fastq    = $working_dir."/".$project_id.".fq";
`cp $fasta_in $ref` unless (-e $ref);

`$bwa index -a is $ref`; 
`$samtools faidx $ref`;

###(2)The reads were mapped to the reference using BWA.
print STDERR "\n\nStart bwa mem\n\n";

my $fastq = $fastq_in;
$fastq =~ s/,/ /g;
$cmd = "$bwa mem $ref $fastq > $working_dir/aln.sam";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.sam");

###(3) convert bam file to sam file
print STDERR "\n\nStart samtool view\n\n";

$cmd = "$samtools view -bS $working_dir/aln.sam > $working_dir/aln.bam";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.bam");
 
###(4) sort the bam file
print STDERR "\n\nStart samtools sort\n\n";
$cmd = "$samtools sort $working_dir/aln.bam $working_dir/aln.sorted";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.sorted");

###(5) for sorted bam file, build index
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

###2. convert bam file to generate bcf file
print STDERR "\n\nStart samtools mpileup\n\n";
$cmd = "$samtools mpileup -guSDf $ref $working_dir/aln.sorted.rmdup.bam > $working_dir/aln.bcf";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/aln.bcf");

###3. call bcftools to call SNP and INDEL
print STDERR "\n\nStart bcftools\n\n";
$cmd = "$bcftools view $working_dir/aln.bcf > $working_dir/snp_indel.vcf";
print STDERR $cmd, "\n\n";
`$cmd` unless (-e "$working_dir/snp_indel.vcf");

print STDERR "\n\nCompleted!\n\n";


#!/usr/bin/perl

use strict;
use warnings;

use FindBin qw( $Bin );
my $project_id     = shift || "1234";
my $working_dir    = shift || "/tmp";
my $fasta_in       = shift || "";
my $fastq_in       = shift || "";
	 

print STDERR "\nperl rnaedit_disc.pl 1234 /tmp/o40 sample_files/sample.fas sample_files/SRR1004791_1_6w.fq,sample_files/SRR1004791_2_6w.fq\n";
print STDERR "$project_id\t$fasta_in\t$fastq_in\n";
die "\nUsage: perl $0 project_id working_dir file_fasta file_fastq1,file_fastq2\n\n" unless ($project_id && -e $fasta_in);
#die "\nUsage: perl $0 working_dir project_id file_fasta file_fastq1,file_fastq2\n\n" unless ($project_id && -e $fasta_in && -e $fastq_in);

#my $working_dir = "/tmp/dir_".$project_id;
$working_dir =~ s/\/+$//;
my $index    = $working_dir."/".$project_id;
`mkdir $working_dir` unless (-e $working_dir);
my $fasta   = $working_dir."/".$project_id.".fa";
my $fastq   = $working_dir."/".$project_id.".fq";
`cp $fasta_in $fasta` unless (-e $fasta);
#`cp $fastq_in $fastq` unless (-e $fastq);

###RNA edit site analysis 
###1.mapping reads to reference genome

#my $python_cmd = "/home/cliu/anaconda2/bin/python2.7";
my $python_cmd = "/home/cliu/anaconda2/envs/rnaedit/bin/python";
my $samtools   = $Bin."/exe/samtools-0.1.18/samtools";
my $tophat     = $Bin."/exe/tophat-2.1.1.Linux_x86_64/tophat";
#my $bt_HOME    = $Bin."/../bowtie-0.12.7";
#my $bt_HOME   = "/home/cliu/workdir/my_apps/bowtie-0.12.7";
my $bt2_HOME   = $Bin."/exe/bowtie2-2.3.4.3/";
#my $bt2_HOME = "/home/cliu/my_apps/bowtie2-2.3.4.3";
print STDERR "\nrunning Bowtie2\n\n";

my $tophat_out = $working_dir."/dir_th.out";
`mkdir $tophat_out` unless (-e $tophat_out);

my $cmd = "$bt2_HOME/bowtie2-build $fasta $index";
#my $cmd = "$bt_HOME/bowtie-build $fasta $index";
print STDERR "\n$cmd\n\n";
`$cmd`;

print STDERR "\nrunning tophat\n\n";
#$cmd = "$tophat -p 8 -o $tophat_out $index $fastq";
$fastq = $fastq_in;
$fastq =~ s/,/ /g;
$cmd = "(export PATH=$bt2_HOME:\$PATH && $tophat -p 8 -o $tophat_out --library-type=fr-firststrand --read-edit-dist 7 $index $fastq)";
print STDERR "\n$cmd\n\n";
`$cmd` unless (-e $tophat_out."/accepted_hits.bam");
#`tophat -p 8 -o $tophat_out --library-type=fr-firststrand --read-edit-dist 7 $index $fastq` unless (-e $tophat_out."/accepted_hits.bam");

my $bam_forward_sorted = $tophat_out."/f_sorted";
my $bam_reverse_sorted = $tophat_out."/r_sorted";

###2分别得到map到postive strand 和negative strand 上的reads

print STDERR "Extract reads mapped to the positive strand and negative strand\n";

$cmd = "$samtools faidx $fasta";
print STDERR "\n$cmd\n\n";
`$cmd`;

$cmd = "$samtools view -F 0x10 $tophat_out/accepted_hits.bam | $samtools view -bS -T $fasta - | $samtools sort - $bam_forward_sorted";
print STDERR "\n$cmd\n\n";
`$cmd` unless (-e $bam_forward_sorted);

$cmd = "$samtools view -f 0x10 $tophat_out/accepted_hits.bam | $samtools view -bS -T $fasta - | $samtools sort - $bam_reverse_sorted";
print STDERR "\n$cmd\n\n";
`$cmd` unless (-e $bam_reverse_sorted);

#`samtools view $tophat_out/accepted_hits.bam | gawk '(! and(16, $2))' | samtools vie2yyw -bS -T $fasta - | samtools sort - $bam_forward_sorted`;
#`samtools view $tophat_out/accepted_hits.bam | gawk '(and(16, $2))' | samtools view -bS -T $fasta - | samtools sort - $bam_reverse_sorted`;
###3用REDItoolDenovo.py 预测RNA编辑位点。
$bam_forward_sorted = $bam_forward_sorted.".bam";
$bam_reverse_sorted = $bam_reverse_sorted.".bam";
my $out_f = $working_dir."/forward_out.txt";
my $out_r = $working_dir."/reverse_out.txt";

$cmd = "$python_cmd $Bin/exe/REDItools-1.0.4/reditools/REDItoolDenovo.py -i $bam_forward_sorted -f $fasta -s 12 -c 5 -l -o $out_f -t 7 -v 1";
print STDERR "\n$cmd\n\n";
`$cmd`;

$cmd = "$python_cmd $Bin/exe/REDItools-1.0.4/reditools/REDItoolDenovo.py -i $bam_reverse_sorted -f $fasta -s 12 -c 5 -l -o $out_r -t 7 -v 1";
print STDERR "\n$cmd\n\n";
`$cmd`; 

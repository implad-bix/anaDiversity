#encoding=utf-8
import sys
from Bio import SeqIO
import os
import subprocess
import fnmatch
import getopt

homepath = os.path.dirname(os.path.abspath(sys.argv[0]))
blast_path = os.path.abspath(homepath+"/../exe")

################################################
#print >>sys.stderr, "version: v20190412"
#print >>sys.stderr, "usage: queryFile queryType[fasta|gb] leftFastqFile rightFastqFile outPrefix max_target_seqs outTmpDir"
#print >>sys.stderr, "example: python2.7 filter_reads_for_plastome_RNAEditingAnalysis.py NC_000932.gb gb reads_1.fq reads_2.fq select 100 outTmpDir"
#print "-----------------------------------------------------------------------\n"
#opts, args = getopt.getopt(sys.argv[1:], "hd:q:o:B", ["help", "database=", "query=", "outdir=", "Bootstrap"])

opts, args = getopt.getopt(sys.argv[1:], "hvq:t:l:r:p:m:o:", ["help","version", "query=", "type=", "left=", "right=", "prefix=", "max_depth=100", "outdir="])



def usage():
	print """\nPREA 0.1 (April 2019)
Copyright (C) 2019 Linchun Shi and Chang Liu
      		Institute of Medicinal Plant Development
  		Chinese Academy of Medical Science
Freely distributed under the GNU General Public License (GPLv3)
"""
	print "Usage: python2.7 -q queryFile -t queryType[fasta|gb] -l leftFastqFile -r rightFastqFile -p outPrefix -m max_depth -o outdir\n"
	print "Example: python2.7 -q NC_000932.gb -t gb -l reads_1.fq -r reads_2.fq -p Arth -m 100 -o outTmpDir\n"
	print "-h, --help	print help"
	print "-v, --version	version"
	print "-q, --query 	query file"
	print "-t, --type	query file type, it should be gb or fasta"
	print "-l, --left	the left reads file"
	print "-r, --right	the right reads file"
	print "-p, --prefix	the prefix name for output file"
	print "-m, --max_depth	the max depth for genes"
	print "-o, --outdir	the outdir for depositing temporay files\n"
	print """
The output format will be prefix_top_max_depth_1/2.fq.\n
"""
queryFile = ""
queryType = ""
leftFastqFile = ""
rightFastqFile = ""
outPrefix = ""
max_target_seqs = "100"
outTmpDir = ""

if not  opts:
	usage()
	sys.exit()
else:
    for o, a in opts:
	if o in ("-h", "--help"):
		usage()
		sys.exit()
	if o in ("-v", "--version"):
		print "[*] Version is v20190412"
		sys.exit()
	if o in ("-q", "--query"):
		queryFile = a
	if o in ("-t", "--type"):
		queryType = a
	#if o in ("-", "--"):
	if o in ("-l", "--left"):
		leftFastqFile = a
	if o in ("-r", "--right"):
		rightFastqFile = a
	if o in ("-p", "--prefix"):	
		outPrefix = a
	if o in ("-m", "--max_depth"):
		max_target_seqs = a
	if o in ("-o", "--outdir"):
		outTmpDir = a

if not queryFile:
	print "please provide queryFile"
	sys.exit()
if not queryType:
	print "please specify queryType"
	sys.exit()
if not leftFastqFile:
	print "please provide left reads in fastq format"
	sys.exit()
if not rightFastqFile:
	print "please provide right reads in fastq format"
        sys.exit()
if not outPrefix:
	print "please provide prefix name"
	sys.exit()
if not outTmpDir:
	print "pleast specify output dir"
	sys.exit()

#queryFile = sys.argv[1] #-q --query 
#queryType = sys.argv[2] #fasta|gb -t --type 
#leftFastqFile = sys.argv[3] #-1 --left
#rightFastqFile = sys.argv[4] #-2 --right
#outPrefix = sys.argv[5] #-p --prefix
#max_target_seqs = sys.argv[6] #-m --max_target_seqs
#outTmpDir = sys.argv[7] #-o --outdir

forwardFastqFile = leftFastqFile
reverseFastqFile = rightFastqFile

outHandle = outPrefix
outdir = outTmpDir

if not os.path.exists(outdir):
	os.mkdir(outdir)

#################################################
def extractGene8CDSSequence(gbFile, outfile):
    out = open(outfile, "w")
    L_CDS_gene = []
    workfile = gbFile

    for record in SeqIO.parse(workfile, "gb"):
        for feature in record.features:
	    if feature.type == "CDS":
		if feature.qualifiers.has_key('gene'):
			cds_gene_name = feature.qualifiers['gene'][0]
		else:
			cds_gene_name = "ngn" # no gene name
		L_CDS_gene.append(cds_gene_name)

    #print L_CDS_gene

    for record in SeqIO.parse(workfile, "gb"):
	m_gene = 0
        m_cds = 0
	accession = record.annotations['accessions'][0].replace("_","-")
	organism = record.annotations['organism'].replace(" ", "_")
	flag = organism.split("_")[0][0:2]+organism.split("_")[1][0:2]+"Cp"+accession

        for feature in record.features:
	    if feature.type == "gene":
		if feature.qualifiers.has_key('gene'):
			gene_name = feature.qualifiers['gene'][0]
		else:
			gene_name = "ngn"

		if gene_name in L_CDS_gene:
			m_gene += 1
			gene_locus_tag = flag+"_gene"+str(m_gene)
			gene_sequence = str(feature.extract(record.seq))
			out.write(">"+gene_locus_tag+"_"+gene_name+"\n"+gene_sequence+"\n")

	    if feature.type == "CDS":
		if feature.qualifiers.has_key('gene'):
                        cds_gene_name = feature.qualifiers['gene'][0]
                else:
                        cds_gene_name = "ngn" # no gene name
		m_cds += 1
		cds_locus_tag = flag+"_cds"+str(m_cds)
		cds_sequence = str(feature.extract(record.seq))
		out.write(">"+cds_locus_tag+"_"+cds_gene_name+"\n"+cds_sequence+"\n")
    out.close()

outQueryFastaFile = outdir+"/query.fasta"

if queryType == "fasta":
	L_record = []
	for record in SeqIO.parse(queryFile, "fasta"):
		L_record.append(record)

	if L_record:
		SeqIO.write(L_record, outQueryFastaFile, "fasta")
	else:
		print >>sys.stderr, "No fasta record"
		sys.exit()

elif queryType == "gb":
	L_record = []
	for record in SeqIO.parse(queryFile, "gb"):
		L_record.append(record)
	if L_record:
		extractGene8CDSSequence(queryFile, outQueryFastaFile)

	else:
		print >>sys.stderr, "No gb record"
		sys.exit()
else:
	print >>sys.stderr, "queryType error, it should be word \"fasta\" or \"gb\""
	sys.exit()

#outHandle = "blastdb"
forwardFastaFile = outdir+"/"+outHandle+"_1.fasta"
reverseFastaFile = outdir+"/"+outHandle+"_2.fasta"

forwardBlastoutFile = outdir+"/"+outHandle+"_1.B.O"
reverseBlastoutFile = outdir+"/"+outHandle+"_2.B.O"

outSbjctFile = outdir+"/"+outHandle+"_sbjct.txt"

outFowradTargetFastqFile = outdir+"/"+outHandle+"_top"+max_target_seqs+"_1.fq"
outReverseTargetFastqFile = outdir+"/"+outHandle+"_top"+max_target_seqs+"_2.fq"
if not os.path.exists(forwardFastaFile):
	cmd1 = "%s/seqkit fq2fa -j 8 %s -o %s" % (homepath+"/../exe", forwardFastqFile, forwardFastaFile)
	P = subprocess.Popen(cmd1, shell=1)
	P.wait()
else:
	print "+++\tSkip convert leftFastq to leftFasta\n"

#reverseFastaFile = outdir+"/"+outHandle+"_2.fasta"
if not os.path.exists(reverseFastaFile):
	cmd2 = "%s/seqkit fq2fa -j 8 %s -o %s" % (homepath+"/../exe",reverseFastqFile, reverseFastaFile)
	P = subprocess.Popen(cmd2, shell=1)
	P.wait()
else:
	print "+++\tSkip convert rightFastq to rightFasta\n"

import fnmatch

if not fnmatch.filter(os.listdir(outdir), outHandle+"_1.fasta"+"*.nin"):
	cmd3 = "%s/makeblastdb -in %s -dbtype nucl -parse_seqids" % (blast_path, forwardFastaFile)
	P = subprocess.Popen(cmd3, shell=1)
	P.wait()
else:
	print "+++\tSkip makeblastdb for leftFastaFile\n"

if not fnmatch.filter(os.listdir(outdir), outHandle+"_2.fasta"+"*.nin"):
	cmd4= "%s/makeblastdb -in %s -dbtype nucl -parse_seqids" % (blast_path, reverseFastaFile)
	P = subprocess.Popen(cmd4, shell=1)
	P.wait()
else:
	print "+++\tSkip makeblastdb for print rightFastaFile\n"

#forwardBlastoutFile = outdir+"/"+outHandle+"_1.B.O"
#reverseBlastoutFile = outdir+"/"+outHandle+"_2.B.O"
 
cmd5 = "%s/blastn -query %s -db %s -outfmt 6 -max_target_seqs %s -evalue 1e-6 -num_threads 8 -task megablast -out %s" % (blast_path, outQueryFastaFile, forwardFastaFile, max_target_seqs, forwardBlastoutFile)
P = subprocess.Popen(cmd5, shell=1)
P.wait()

cmd6 = "%s/blastn -query %s -db %s -outfmt 6 -max_target_seqs %s -evalue 1e-6 -num_threads 8 -task megablast -out %s" % (blast_path, outQueryFastaFile, reverseFastaFile, max_target_seqs, reverseBlastoutFile)
P = subprocess.Popen(cmd6, shell=1)
P.wait()

L_sbjct = []

for line in open(forwardBlastoutFile):
	line = line.strip()
	sbjct = line.split("\t")[1]
	if sbjct not in L_sbjct:
		L_sbjct.append(sbjct)

for line in open(reverseBlastoutFile):
	line = line.strip()
        sbjct = line.split("\t")[1]
        if sbjct not in L_sbjct:
                L_sbjct.append(sbjct)
#outSbjctFile = outdir+"/"+outHandle+"_sbjct.txt"
outSbjct = open(outSbjctFile, "w")

for sbjct in L_sbjct:
	outSbjct.write(sbjct+"\n")
outSbjct.close()

#outFowradTargetFastqFile = outHandle+"_top"+max_target_seqs+"_1.fq"
#outReverseTargetFastqFile = outHandle+"_top"+max_target_seqs+"_2.fq"

cmd_notCombined_1 = "%s/seqtk subseq %s %s >%s" % (homepath+"/../exe/seqtk-1.3", forwardFastqFile, outSbjctFile, outFowradTargetFastqFile)
P = subprocess.Popen(cmd_notCombined_1, shell=1)
P.wait() 

cmd_notCombined_2 = "%s/seqtk subseq %s %s >%s" % (homepath+"/../exe/seqtk-1.3", reverseFastqFile, outSbjctFile, outReverseTargetFastqFile)
P = subprocess.Popen(cmd_notCombined_2, shell=1)
P.wait()


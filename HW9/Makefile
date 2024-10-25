GENOME=GCF_910591885.1
SRA=SRR3693403
N=1000000

REF=ref/b_terrestris.fna
R1=reads/${SRA}_1.fastq
R2=reads/${SRA}_2.fastq

BAM=bam/b_terrestris.bam
Q10=bam/q10/b_terrestris_q10.bam

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
================================================================================================
info: 
	@echo "make download # download the references genome"
	@echo "make index # index the reference genome"
	@echo "make align # align the reads"
	@echo "make flagstat # generate the diferent flag stat for the bam file"
	@echo "make q10 # create a new bam file with a minimum mapping quality of 10"
	@echo "make clean # remove the files"

download:
	# make the directory for the reference genome
	mkdir -p $(dir ${REF})
	# download the reference genome
	datasets download genome accession ${genome}
	# unzip the file
	unzip -n ncbi_dataset.zip -x README.md md5sum.txt
	# copy the reference genome to the ref directory
	cp -f ncbi_dataset/data/${genome}/*_genomic.fna ${REF}
	# remove the ncbi_dataset directory
	rm -rf ncbi_dataset

fastq:
	# make the directory for the reads
	mkdir -p $(dir ${R1})
	# download the SRA file
	fastq-dump -X ${N} --outdir reads --split-files ${SRA} 

index: 
	mkdir -p ref
	bwa index ${REF}

align: 
	# make the directory for the bam file
	mkdir -p $(dir ${BAM})
	# align the reads
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort > ${BAM}
	# index the bam file
	samtools index ${BAM}

flagstat:
	# number of reads that don't aligned to the reference
	samtools view -c -f 4 ${BAM}
	# number of primary reads aligned to the reference
	samtools view -c -F 256 ${BAM}
	# number of secondary reads aligned to the reference
	samtools view -c -f 256 ${BAM}
	# number of supplementary reads aligned to the reference
	samtools view -c -f 2048 ${BAM}
	# number of properly-paired alignments on the reverse strand
	samtools view -c -f 99 ${BAM}

q10: 
	#make the directory for the q10 bam file
	mkdir -p $(dir ${Q10})
	# create a new bam file with a minimum mapping quality of 10
	samtools view -h -b -q 10 -f 1 -F 256 ${BAM} > ${Q10}
	samtools index bam/b_terrestris_q10.bam
	# compare the flagstat of the original bam file and the q10 bam file
	samtools flagstat ${BAM}
	samtools flagstat ${Q10}

clean:
	rm -rf ref
	rm -rf reads
	rm -rf bam
	rm -rf ncbi_dataset.zip
	rm -rf ncbi_dataset

.PHONY: download fastq index align flagstat q10 

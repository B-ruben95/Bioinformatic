
# Accession number of Bombus terrestris genome.
GENOME=GCF_910591885.1
# The sequencing read accession number.
SRA?=
# The number of reads to get
N=1000000
# The reference file
REF=ref/b_terrestris.fna
# The path to read 1
R1=reads/${SRA}_1.fastq
# The path to read 2
R2=reads/${SRA}_2.fastq
# The GFF file.
GFF=ref/b_terrestris.gff
# The name of the sample
SAMPLE?=I-D1
# The resulting BAM file.
BAM=bam/b_terrestris_${SRA}.bam
# The resulting variant VCF file
VCF=vcf/b_terrestris_${SRA}.vcf.gz

# Additional flags for mpileup.
PILE_FLAGS =  -d 100 --annotate 'INFO/AD,FORMAT/DP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/SP'
# Additional bcf flags for calling process.
CALL_FLAGS = --ploidy 2 --annotate 'FORMAT/GQ'

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules
###########################################################################################################################################
info: 
	@echo "make download # download the references genome"
	@echo "make fastq # download the reads"
	@echo "make index # index the reference genome"
	@echo "make align # align the reads"
	@echo "make flagstat # generate the different flag stat for the bam file"
	@echo "make q10 # create a new bam file with a minimum mapping quality of 10"
	@echo "make clean # remove the files"

check_files:
	@if [ ! -f ref/b_terrestris.fna ] || [ ! -f ref/b_terrestris.fna.bwt ]; then \
		echo "#"; \
		echo "# Please download and index the genome first with: make download index"; \
		echo "#"; \
		exit 1; \
	fi

all: check_files fastq align variant


all: fastq align variant

download:
	# make the directory for the reference genome
	mkdir -p $(dir ${REF})
	# download the reference genome 
	datasets download genome accession ${GENOME} --include gff3,genome
	# unzip the file
	unzip -n ncbi_dataset.zip -x README.md md5sum.txt
	# copy the reference genome to the ref directory
	cp -f ncbi_dataset/data/${GENOME}/*_genomic.fna ${REF}
	# copy the gff file to the ref directory
	cp -f ncbi_dataset/data/${GENOME}/genomic.gff ${GFF}
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

variant:
	# make the directory for the vcf file
	mkdir -p $(dir ${VCF})
	# call the variants
	bcftools mpileup ${PILE_FLAGS} -O u -f ${REF} ${BAM} | \
		bcftools call ${CALL_FLAGS} -mv -O u | \
		bcftools norm -f ${REF} -d all -O u | \
		bcftools sort -O z > ${VCF}
	# index the vcf file
	bcftools index ${VCF}
	# generate the stats report
	bcftools stats ${VCF} > stats_report.txt

merge:
	# Merge VCF files into a single one.
	bcftools merge -0 vcf/*.vcf.gz -O z > merged.vcf.gz

	# Index the merged VCF file
	bcftools index merged.vcf.gz

stats:	
	@grep "SN[[:space:]]0[[:space:]]number of samples" merged_stats_report.txt | awk '{print "    Samples:", $$6}'
	@grep "SN[[:space:]]0[[:space:]]number of SNPs" merged_stats_report.txt | awk '{print "    SNPs:", $$6}'
	@grep "SN[[:space:]]0[[:space:]]number of indels" merged_stats_report.txt | awk '{print "    Indels:", $$6}'
	@grep "SN[[:space:]]0[[:space:]]number of records" merged_stats_report.txt | awk '{print "    Total Variants:", $$6}'
	@gunzip -c merged.vcf.gz | grep -v "^#" | awk '{ if ($$6 > 30) print }' | wc -l | awk '{print "    Variants > 30:", $$1}'


clean:
	rm -rf ref
	rm -rf reads
	rm -rf bam
	rm -rf vcf
	rm -rf ncbi_dataset.zip
	rm -rf ncbi_dataset

.PHONY: download fastq index align variant vep

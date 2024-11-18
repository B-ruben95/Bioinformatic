GENOME=GCF_910591885.1
SRA=SRR3693403
N=1000000

REF=ref/b_terrestris.fna
R1=reads/${SRA}_1.fastq
R2=reads/${SRA}_2.fastq

BAM=bam/b_terrestris.bam

VCF=vcf/b_terrestris.vcf.gz
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
	@echo "make index # index the reference genome"
	@echo "make align # align the reads"
	@echo "make flagstat # generate the diferent flag stat for the bam file"
	@echo "make q10 # create a new bam file with a minimum mapping quality of 10"
	@echo "make clean # remove the files"

download:
	# make the directory for the reference genome
	mkdir -p $(dir ${REF})
	# download the reference genome
	datasets download genome accession ${GENOME}
	# unzip the file
	unzip -n ncbi_dataset.zip -x README.md md5sum.txt
	# copy the reference genome to the ref directory
	cp -f ncbi_dataset/data/${GENOME}/*_genomic.fna ${REF}
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
	echo "SNPs: $(grep "SN[[:space:]]0[[:space:]]number of SNPs" stats_report.txt | awk '{print $6}')"
	echo "Indels: $(grep "SN[[:space:]]0[[:space:]]number of indels" stats_report.txt | awk '{print $6}')"
	echo "Total Variants: $(grep "SN[[:space:]]0[[:space:]]number of records" stats_report.txt | awk '{print $6}')"
	echo "variants with a quality score greater than 30: $(gunzip -c ${VCF} | grep -v "^#" | awk '$6 > 30' | wc -l)" 
	#variant per chromosome
	gunzip -c ${VCF}| grep -v "^#" | cut -f1 | sort | uniq -c
	#SNP per chromosome
	gunzip -c ${VCF} |grep -v "^#" | awk 'length($4) == 1 && length($5) == 1' | cut -f1 | sort | uniq -c
	#Indels per chromosome
	gunzip -c ${VCF} | grep -v "^#"| awk 'length($4) != 1 || length($5) != 1' | cut -f1 | sort | uniq -c

clean:
	rm -rf ref
	rm -rf reads
	rm -rf bam
	rm -rf vcf
	rm -rf ncbi_dataset.zip
	rm -rf ncbi_dataset

.PHONY: download fastq index align variant 
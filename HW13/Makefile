# Accession number of Bombus terrestris genome.
GENOME=GCF_910591885.1
#the bioproject number
PRJ?=PRJNA1036000
# desing file
DESING=design.csv
# The sequencing read accession number.
SRA?=SRR26724898
# The number of reads to get
N=1000000
# The reference file
REF=ref/b_terrestris.fna
# The path to read 1
R1=reads/${SRA}_1.fastq
# The path to read 2
R2=reads/${SRA}_2.fastq
# The GTF file.
GTF=ref/b_terrestris.gtf
# The name of the sample
SAMPLE?=I-D1
# The resulting BAM file.
BAM=bam/b_terrestris_${SAMPLE}.bam
# The resulting variant VCF file
VCF=vcf/b_terrestris_${SRA}.vcf.gz
# The counts in tab delimited format.
COUNTS_TXT = res/counts-hisat.txt
# Final combinted counts in CSV format.
COUNTS = res/counts-hisat.csv

SHELL = bash
.ONESHELL:
.SHELLFLAGS = -eu -o pipefail -c
.DELETE_ON_ERROR:
MAKEFLAGS += --warn-undefined-variables
MAKEFLAGS += --no-builtin-rules

###########################################################################################################################################

info: 
	@echo "make genome # download the references genome"
	@echo "make fastq # download the reads"
	@echo "make index # index the reference genome"
	@echo "make align # align the reads"
	@echo "make flagstat # generate the different flag stat for the bam file"
	@echo "make q10 # create a new bam file with a minimum mapping quality of 10"
	@echo "make clean # remove the files"

genome:
	# make the directory for the reference genome
	mkdir -p $(dir ${REF})
	# download the reference genome 
	datasets download genome accession ${GENOME} --include gtf,genome
	# unzip the file
	unzip -n ncbi_dataset.zip -x README.md md5sum.txt
	# copy the reference genome to the ref directory
	cp -f ncbi_dataset/data/${GENOME}/*_genomic.fna ${REF}
	# copy the gff file to the ref directory
	cp -f ncbi_dataset/data/${GENOME}/genomic.gtf ${GTF}
	# remove the ncbi_dataset directory
	rm -rf ncbi_dataset

index: 
	mkdir -p ref
	hisat2-build --threads ${NCPU} ${REF} ref/b_terrestris_index

design:
	# Create the design file
	bio search ${PRJ} -H --csv > design.csv
	
rnaseq: fastq align

fastq:
	# make the directory for the reads
	mkdir -p $(dir ${R1})
	# download the SRA file
	fastq-dump -X ${N} --outdir reads --split-files ${SRA} 

align: 
	mkdir -p $(dir ${BAM})
	if [ -z "${R2}" ] || [ ! -s "${R2}" ]; then
		echo "Using single-end mode..."
		hisat2 -p 4 -x ref/b_terrestris_index -U ${R1} | samtools sort -o ${BAM}
	else
		echo "Using paired-end mode..."
		hisat2 -p 4 -x ref/b_terrestris_index -1 ${R1} -2 ${R2} | samtools sort -o ${BAM}
	fi
	samtools index ${BAM}
	bedtools genomecov -ibam ${BAM} -bg > ${BAM}.bedgraph
	samtools faidx ${REF}
	bedGraphToBigWig ${BAM}.bedgraph ${REF}.fai ${BAM}.bw
	rm ${BAM}.bedgraph

count: ${COUNTS_TXT} 

${COUNTS_TXT}:
	# Make the directory name for the counts
	mkdir -p $(dir $@)

	# Count the features
	ls bam/*.bam | xargs featureCounts -a ${GTF} -o $@

csv: ${COUNTS}

${COUNTS}: ${COUNTS_TXT}
	micromamba run -n stats Rscript src/r/format_featurecounts.r -c ${COUNTS_TXT} -o ${COUNTS}

clean:
	rm -rf ref reads bam res vcf design.csv


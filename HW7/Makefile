
#Variables

#variables of the homework 5
NCBI=GCF_910591885.1
TYPE=genome
GENOME=ncbi_dataset/data/GCF_910591885.1/GCF_910591885.1_iyBomTerr1.2_genomic.fna
CHROMOSOME7=chromosome7.fa
R1=reads/wgsim_read1.fq
R2=reads/wgsim_read2.fq

#variables of the homework 6
SRA=SRR12432506
N=100000
FASTQ1=READS/${SRA}_1.fastq
FASTQ2=READS/${SRA}_2.fastq

TRIMMED1=READS/${SRA}_1.fastq.trimmed.fastq
TRIMMED2=READS/${SRA}_2.fastq.trimmed.fastq

usage:
	@echo "make genome  # download the genome file"
	@echo "make stadistics # stadistics of the genome"
	@echo "make simulate  # simulate reads from chromosome 7"
	@echo "make download  # download the data from the SRA database"
	@echo "make quality   # quality control"
	@echo "make trim      # trimming"
	@echo "make multiqc   # multiqc"
	@echo "make clean     # remove the downloaded files"

#Download the genome
genome:
	ncbi-genome-download --annotation-type ${TYPE} ${NCBI} --include ${TYPE}
	@echo "The genome has been downloaded"
	unzip ncbi_dataset.zip

stadistics:
	SIZE=$(du -h ${GENOME} | cut -f1)
	@echo "The size of the genome is ${SIZE}"
	TOTAL_SIZE=$(grep -v ">" ${GENOME} | wc -c)
	@echo "The total size of the genome is ${TOTAL_SIZE}"
	CHROMOSOMES=$(grep -c ">" ${GENOME})
	@echo "The number of chromosomes is ${CHROMOSOMES}"
	awk '/^>/ {if (seqlen){print seqname, seqlen}; seqname=$1; sub(/^>/, "", seqname); seqlen=0; next} {seqlen += length($0)} END {print seqname, seqlen}' ${GENOME} 

simulate:
	cat ${GENOME} | sed -n '/^>NC_063275.1/,/^>/p' | sed '$d' > chromosome7.fa
	mkdir -p $(dirname ${R1})
	wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 ${CHROMOSOME7} ${R1} ${R2}
	seqkit stats ${R1} ${R2}
	READS=$(grep -c "@" ${R1})
	@echo "The number of reads is ${READS}"
	AVG_LENGTH=$(seqkit stats ${R1} ${R2} | grep -v "file" | awk '{s+=$6} END {print s/NR}')
	@echo "The average length of the reads is ${AVG_LENGTH}"
	SIZE_R1=$(du -h ${R1} | cut -f1)
	SIZE_R2=$(du -h ${R2} | cut -f1)
	@echo "The size of the reads is ${SIZE_R1} and ${SIZE_R2}"
 
directory:
	@echo make a directory for the data
	mkdir -p READS REPORTS

download:
	#Download the data
	fastq-dump -X ${N} --split-files ${SRA} -O READS
	@echo "The data have been downloaded"

quality:
	#Quality control
	fastqc -q -o REPORTS ${FASTQ1} ${FASTQ2}
	@echo "Quality control done"

trim:
	fastp -i ${FASTQ1} -I ${FASTQ2} -o ${TRIMMED1} -O ${TRIMMED2} \
      --detect_adapter_for_pe --cut_tail -w 4
	@echo "Trimming done"
	#Quality control
	fastqc -q -o REPORTS ${TRIMMED1} ${TRIMMED2}
	@echo "Quality control done"

multiqc:
	micromamba run -n menv multiqc -o READS REPORTS 
	@echo "Multiqc done"

clean:
	rm -rf READS REPORTS
	@echo "All files have been removed"

.PHONY: genome stadistics simulate download quality trim multiqc 












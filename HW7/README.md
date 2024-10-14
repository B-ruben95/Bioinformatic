# Week 7: Write a Makefile
**Ruben Martin 

In this part of the file i define the variable that will be use, I divide the variable in the diferent homeworks
```bash
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
```

Here i discribe the difirent targest that will run in and what they do 
``` bash
usage:
	@echo "make genome  # download the genome file"
	@echo "make stadistics # stadistics of the genome"
	@echo "make simulate  # simulate reads from chromosome 7"
	@echo "make download  # download the data from the SRA database"
	@echo "make quality   # quality control"
	@echo "make trim      # trimming"
	@echo "make multiqc   # multiqc"
	@echo "make clean     # remove the downloaded files"
```

This target will downaload the genome 
```bash
genome:
	datasets download genome accession ${NCBI} --include ${TYPE}
	@echo "The genome has been downloaded"
	unzip ncbi_dataset.zip
```

This target will generate some stadistics about the genome 
```bash
stadistics:
	SIZE=$(du -h ${GENOME} | cut -f1)
	@echo "The size of the genome is ${SIZE}"
	TOTAL_SIZE=$(grep -v ">" ${GENOME} | wc -c)
	@echo "The total size of the genome is ${TOTAL_SIZE}"
	CHROMOSOMES=$(grep -c ">" ${GENOME})
	@echo "The number of chromosomes is ${CHROMOSOMES}"
	awk '/^>/ {if (seqlen){print seqname, seqlen}; seqname=$1; sub(/^>/, "", seqname); seqlen=0; next} {seqlen += length($0)} END {print seqname, seqlen}' ${GENOME} 
```

This target will generate the simulation for one of the chromosomes of the genome

```bash
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
```

This target creata two directories to store the data and download the data for the SRA database
```bash
download:
	@echo make a directory for the data
	mkdir -p READS REPORTS
    #Download the data
	fastq-dump -X ${N} --split-files ${SRA} -O READS
	@echo "The data have been downloaded"
```

This target will generate the quality control for the data
```bash
quality:
	#Quality control
	fastqc -q -o REPORTS ${FASTQ1} ${FASTQ2}
	@echo "Quality control done"
```

This target will cut the adapters and generate reads of best quality 
```bash
trim:
	fastp -i ${FASTQ1} -I ${FASTQ2} -o ${TRIMMED1} -O ${TRIMMED2} \
      --detect_adapter_for_pe --cut_tail -w 4
	@echo "Trimming done"
	#Quality control
	fastqc -q -o REPORTS ${TRIMMED1} ${TRIMMED2}
	@echo "Quality control done"
```
This target will generate a multiqc analysis 
```bash
multiqc:
	micromamba run -n menv multiqc -o READS REPORTS 
	@echo "Multiqc done"
```
This target will clean the directory 
```bash
clean:
	rm -rf READS REPORTS ncbi_dataset
	@echo "All files have been removed"
```
This target will exectute the hold script

```bash
.PHONY: genome stadistics simulate download quality trim multiqc
```


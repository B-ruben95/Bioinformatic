#variables

genome=GCF_910591885.1
SRA=SRR3693403

REF=ref/b_terrestris.fna
CH=ch/b_terrestris_ch7.fna
R1=reads/${SRA}_1.fastq
R2=reads/${SRA}_2.fastq
RS1=reads/simulated/b_terrestris_ch7.fastq
RS2=reads/simulated/b_terrestris_ch7_2.fastq
N=10000000
BAM=bam/b_terrestris.bam
BAMS=bam/simulate/b_terrestris_ch7.bam
info:
	@echo

download:

	#make the directory for the reference genome
	#mkdir -p $(dir ${REF})
	#make the directory for the reads
	mkdir -p $(dir ${R1})
	#download the reference genome
	#datasets download genome accession ${genome}
	#unzip the file
	#unzip -n ncbi_dataset.zip -x README.md md5sum.txt
	#download the SRA file
	fastq-dump -X ${N}
	3 --outdir reads --split-files ${SRA}  #i deleated the -X ${N} because the aligment was to smal 
	#copy the reference genome to the ref directory
	#cp -f ncbi_dataset/data/${genome}/*_genomic.fna ${REF}
	
index: 
	mkdir -p ref
	bwa index ${REF}

align: 
	#make the directory for the bam file
	mkdir -p $(dir ${BAM})
	#align the reads
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort > ${BAM}
	#index the bam file
	samtools index ${BAM}

simulate:
	#make the directory for the chromosome 7
	mkdir -p $(dir ${CH})
	#make the directory for the simulated reads
	mkdir -p reads/simulated
	#extract the chromosome 7 from the reference genome
	samtools faidx ${REF} NC_063275.1 > ${CH}
	#simulate the reads
	wgsim -N ${N} -1 100 -2 100 -r 0 -R 0 -X 0 ${CH} ${RS1} ${RS2}
	#make the directory for the bam file
	mkdir -p $(dir ${BAMS})
	#align the reads
	bwa mem -t 4 ${REF} ${RS1} ${RS2} | samtools sort > ${BAMS}
	#index the bam file
	samtools index ${BAMS}

stats: 
	#generate the stats for the bam file
	samtools flagstat ${BAM}
	#generate the stats for the simulated bam file
	samtools flagstat ${BAMS}
	
clean:
	rm -rf ${REF} ${R1} ${R2} ${BAM} ${BAM}.bai
	
.PHONY: download align stats

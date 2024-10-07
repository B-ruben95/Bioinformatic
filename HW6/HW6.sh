set -uex

#Variables

NCBI=SRR12432506

FASTQ1=READS/${NCBI}_1.fastq
FASTQ2=READS/${NCBI}_2.fastq

TRIMMED1=READS/${NCBI}_1.fastq.trimmed.fastq
TRIMMED2=READS/${NCBI}_2.fastq.trimmed.fastq



#make a directory for the data
mkdir -p READS REPORTS


#Download the data
fastq-dump --split-files $NCBI -O READS
echo "The data have been downloaded"

#Quality control
fastqc -q -o REPORTS $FASTQ1 $FASTQ2
echo "Quality control done"

#Trimming

fastp -i ${FASTQ1} -I ${FASTQ2} -o ${TRIMMED1} -O ${TRIMMED2} \
      --detect_adapter_for_pe --cut_tail -w 4
echo "Trimming done"

#Quality control

fastqc -q -o REPORTS $TRIMMED1 $TRIMMED2
echo "Quality control done"















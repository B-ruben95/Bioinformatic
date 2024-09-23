set -uex

#variables

NCBI=GCF_000005845.2
TYPE=gff3,genome
GENOME=ncbi_dataset/data/${NCBI}/${NCBI}*.fna
GENOMIC=ncbi_dataset/data/$NCBI/genomic.gff

#
# don't change anything below this line
#

# Download the file from the NCBI if the file doesn't exist
[ -f ncbi_dataset.zip ] || datasets download genome accession ${NCBI} --include ${TYPE}

# Unzip the file if it is not already unzipped
[ -d ncbi_dataset ] || unzip ncbi_dataset.zip 

#make fasta index
samtools faidx ${GENOME} > ${GENOME}.fai

#print the number of sequences in the fasta file
cat ${GENOME}.fai | head 


# Print statistics of the genome file
seqkit stats $GENOME

#make a new gff3 file with only the gene information
cat "$GENOMIC" | awk '$3 == "gene"' >genes.gff

#print the number of genes
cat genes.gff | wc -l










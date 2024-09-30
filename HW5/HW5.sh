set -uex

#Variables

NCBI=GCF_910591885.1
TYPE=genome
GENOME=ncbi_dataset/data/GCF_910591885.1/GCF_910591885.1_iyBomTerr1.2_genomic.fna



# Download the genome
[ -f ncbi_dataset.zip ] || datasets download genome accession ${NCBI} --include ${TYPE}

# Unzip the file if it is not already unzipped
[ -d ncbi_dataset ] || unzip ncbi_dataset.zip



#size of the file, total size of the genome, number of chromosomes and the id and length of each chromosome
seqkit stats $GENOME




#number of chromosomes
CHROMOSOMES=$(grep -c ">" ${GENOME})

echo "The number of chromosomes is $CHROMOSOMES"

#name of the chromosomes
CHROM_ID=$(grep ">" ${GENOME} | sed 's/>//g' | sed 's/ .*//g' > chrom_names.txt)

echo "The names of the chromosomes are:"
cat chrom_names.txt


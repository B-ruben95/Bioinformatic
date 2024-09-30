set -ue

#Variables

NCBI=GCF_910591885.1
TYPE=genome
GENOME=ncbi_dataset/data/GCF_910591885.1/GCF_910591885.1_iyBomTerr1.2_genomic.fna

# Part 1
echo 'Bombus terrestris genome statistics'
# Download the genome
[ -f ncbi_dataset.zip ] || datasets download genome accession ${NCBI} --include ${TYPE}

# Unzip the file if it is not already unzipped
[ -d ncbi_dataset ] || unzip ncbi_dataset.zip

# the size of the fasta file
SIZE=$(du -h ${GENOME} | cut -f1)
echo "The size of the fasta file is $SIZE"

#total size of the genome with 
TOTAL_SIZE=$(grep -v ">" ${GENOME} | wc -c)
echo "The total size of the genome is $TOTAL_SIZE"

#number of chromosomes
CHROMOSOMES=$(grep -c ">" ${GENOME})

echo "The number of chromosomes is $CHROMOSOMES"

#name of the chromosomes and their length chromosome in the genome 
awk '/^>/ {if (seqlen){print seqname, seqlen}; seqname=$1; sub(/^>/, "", seqname); seqlen=0; next} {seqlen += length($0)} END {print seqname, seqlen}' $GENOME 


# Part 2
echo 'simulating reads'
#Extract the chromosome 7 from the genome 
cat $GENOME | sed -n '/^>NC_063275.1/,/^>/p' | sed '$d' > chromosome7.fa

CHROMOSOME7=chromosome7.fa

G=26000000
C=10
L=100
N=$(echo "$C * $G / $L" | bc)

#echo "The number of reads needed to cover the genome 10x is $N"

R1=reads/wgsim_read1.fq
R2=reads/wgsim_read2.fq

# Make the directory that will hold the reads extracts 
# the directory portion of the file path from the read
mkdir -p $(dirname ${R1})

# Simulate with no errors and no mutations
wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 ${CHROMOSOME7} ${R1} ${R2}

# Run read statistics
seqkit stats ${R1} ${R2}

# number of reads generated
READS=$(grep -c "@" ${R1})
echo "The number of reads generated is $READS"

# average read length
AVG_LENGTH=$(seqkit stats ${R1} ${R2} | grep -v "file" | awk '{s+=$6} END {print s/NR}')
echo "The average read length is $AVG_LENGTH"

#size of each file 
SIZE_R1=$(du -h ${R1} | cut -f1)
SIZE_R2=$(du -h ${R2} | cut -f1)
echo "The size of the read1 file is $SIZE_R1"
echo "The size of the read2 file is $SIZE_R2"

# compress the files
COMPRESS_R1=$(gzip -k ${R1})
COMPRESS_R2=$(gzip -k ${R2})

#size of the compressed files
SIZE_COMPRESS_R1=$(du -h ${R1}.gz | cut -f1)
SIZE_COMPRESS_R2=$(du -h ${R2}.gz | cut -f1)
echo "The size of the compressed read1 file is $SIZE_COMPRESS_R1"
echo "The size of the compressed read2 file is $SIZE_COMPRESS_R2"

# Space saved
SPACE_SAVED=$(echo "($SIZE_R1 + $SIZE_R2) - ($SIZE_COMPRESS_R1 + $SIZE_COMPRESS_R2)" | bc)
echo "The space saved is $SPACE_SAVED"
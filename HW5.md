# Week 5: Simulating FASTQ files

**Ruben Martin** 

For this homework i used the *Bombus terrestris* genome GCF_910591885.1

### Part 1: Select a genome, then download the corresponding FASTA file.

* The size of the file
```bash
SIZE=$(du -h ${GENOME} | cut -f1)
```
```
The size of the fasta file is 392M
```
* The total size of the genome

```bash
TOTAL_SIZE=$(grep -v ">" ${GENOME} | wc -c)
```
```
The total size of the genome is  397874189
```

* The number of chromosomes in the genome

```bash
CHROMOSOMES=$(grep -c ">" ${GENOME})
```

```bash
The number of chromosomes is 249
```

* The name (id) and length of each chromosome in the genome.

```bash
awk '/^>/ {if (seqlen){print seqname, seqlen}; seqname=$1; sub(/^>/, "", seqname); seqlen=0; next} {seqlen += length($0)} END {print seqname, seqlen}' $GENOME
```

```
NC_063269.1 18372659
NC_063270.1 19735706
NC_063271.1 22099476
NC_063272.1 18668021
NC_063273.1 14194883
NC_063274.1 23129593
NC_063275.1 25524254
NC_063276.1 11657431
NC_063277.1 19446242
NC_063278.1 20491647
NC_063279.1 20781255
NC_063280.1 12897853
NC_063281.1 14593050
NC_063282.1 11964351
NC_063283.1 11413514
NC_063284.1 8887048
NC_063285.1 11646943
NC_063286.1 5271557
```

### Part 2: Generate a simulated FASTQ output for a sequencing instrument of your choice.  Set the parameters so that your target coverage is 10x.

* How many reads have you generated?

```bash
N=$(echo "$EX * $G / $L" | bc)
```

print 
``` bash
The number of reads needed to cover the genome 10x is 2600000
```

* What is the average read length?

```bash
AVG_LENGTH=$(seqkit stats ${R1} ${R2} | grep -v "file" | awk '{s+=$6} END {print s/NR}')
```
Print
```
The average read length is 100
```

* How big are the FASTQ files?

```bash
SIZE_R1=$(du -h ${R1} | cut -f1)
SIZE_R2=$(du -h ${R2} | cut -f1)
```

Print
```
The size of the read1 file is 640M
The size of the read2 file is 640M
```

* Compress the files and report how much space that saves.

``` bash
SIZE_COMPRESS_R1=$(du -h ${R1}.gz | cut -f1)
SIZE_COMPRESS_R2=$(du -h ${R2}.gz | cut -f1)
SPACE_SAVED=$(echo "($SIZE_R1 + $SIZE_R2) - ($SIZE_COMPRESS_R1 + $SIZE_COMPRESS_R2)" | bc)
```
Print
```bash
The size of the compressed read1 file is 129M
The size of the compressed read2 file is 129M
The space saved is 10220
```



### Part 3 : How much data would be generated when covering the Yeast,  the Drosophila or the Human genome at 30x?


$$
Number of reads= \frac {Genome size\ x \ Coverage} {Read length}
$$


| Organism      | Genome Size (bp)      | Number of Reads   | FASTQ Size (Before Compression) | FASTQ Size (After Compression)  |
|:-------------:|:--------------------:|:-----------------:|:-------------------------------:|:-------------------------------:|
| Yeast         | 12,100,000           | 3,630,000         | 1.452 GB                        | 726 MB                          |
| Drosophila    | 143,700,000          | 43,110,000        | 17.244 GB                       | 8.622 GB                        |
| Human         | 3,100,000,000        | 930,000,000       | 372 GB                          | 186 GB                          |








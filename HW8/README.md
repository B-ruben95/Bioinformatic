# Week 8: Generate a BAM alignment file
**Ruben Martin**

For this assigment i uses the references genome for the european buff-tailed bumblebee (*Bombus terrestris*) (GCF_910591885.1) and the SRR from the WGS of this species.

Target for the index
```bash
index: 
	mkdir -p ref
	bwa index ${REF}
```

Target for the aligment 
```bash
align: 
	#make the directory for the bam file
	mkdir -p $(dir ${BAM})
	#align the reads
	bwa mem -t 4 ${REF} ${R1} ${R2} | samtools sort > ${BAM}
	#index the bam file
	samtools index ${BAM}
```
Target to generate the the simulate reads and align them 
```bash
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
```

stats of the reads download from the NCBI 
```bash
samtools index bam/b_terrestris.bam
```
```bash
38156820 + 0 in total (QC-passed reads + QC-failed reads)
37966022 + 0 primary
0 + 0 secondary
190798 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
38049635 + 0 mapped (99.72% : N/A)
37858837 + 0 primary mapped (99.72% : N/A)
37966022 + 0 paired in sequencing
18983011 + 0 read1
18983011 + 0 read2
34584800 + 0 properly paired (91.09% : N/A)
37850570 + 0 with itself and mate mapped
8267 + 0 singletons (0.02% : N/A)
762916 + 0 with mate mapped to a different chr
203701 + 0 with mate mapped to a different chr (mapQ>=5)
```

stats of the simulate reads from the seven chromosomes of B. terrestris
```bash
samtools flagstat bam/simulate/b_terrestris_ch7.bam
```
```bash
2000003 + 0 in total (QC-passed reads + QC-failed reads)
2000000 + 0 primary
0 + 0 secondary
3 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
2000003 + 0 mapped (100.00% : N/A)
2000000 + 0 primary mapped (100.00% : N/A)
2000000 + 0 paired in sequencing
1000000 + 0 read1
1000000 + 0 read2
1999720 + 0 properly paired (99.99% : N/A)
2000000 + 0 with itself and mate mapped
0 + 0 singletons (0.00% : N/A)
40 + 0 with mate mapped to a different chr
7 + 0 with mate mapped to a different chr (mapQ>=5)
```

Reads align to the refereces genome visualises in IGV.

![IGV SRA Data](https://github.com/B-ruben95/Bioinformatic/blob/main/HW8/Image/IGV_B_terrestris.png?raw=true)



* Briefly describe the differences between the two datasets.

The stats for the alimegment show that in both cases the read alignment was successful, with a high proportion of reads mapped and pairs correctly aligned. 
The percentage of mapped reads for the simulated data is higher than for the SRR data, this is mainly because the simulated data was derived directly from the reference genome I used for alignment.
The SRA reads are more than the simulated reads, since only 2 million reads were generated for the simulated reads, so all the statistics were much higher for the simulated data because of the lower number of readings. 


| **Metric**                         | **SRA Reads**             | **Simulated Reads**       | **Comment**                                   |
|------------------------------------|-----------------------------------|----------------------------------|-----------------------------------------------|
| **Total Reads**                    | 38,156,820                       | 2,000,003                        | SRR data contains almost 19x more reads.     |
| **Primary Reads**                  | 37,966,022                       | 2,000,000                        | Both datasets mostly consist of primary reads.|
| **Supplementary Reads**            | 190,798                          | 3                                | SRR data has more split alignments.|
| **Mapped Reads**                   | 38,049,635 (99.72%)              | 2,000,003 (100.00%)              | Simulated data has perfect alignment.         |
| **Properly Paired Reads**          | 34,584,800 (91.09%)              | 1,999,720 (99.99%)               | Simulated data has nearly all reads correctly paired. |
| **Singletons (unpaired reads)**    | 8,267 (0.02%)                    | 0 (0.00%)                        | No singletons are found in the simulated data.|
| **Reads Mapped to Different Chromosomes** | 762,916                       | 40                               | SRR data has more inter-chromosomal alignments. |
| **Reads Mapped to Different Chromosomes (mapQ ≥ 5)** | 203,701 | 7                            | Few reliable inter-chromosomal alignments in simulated data. 


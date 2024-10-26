# Week 9: Filter a BAM file
**Ruben Martin**

For this assigment i uses the references genome for the european buff-tailed bumblebee (*Bombus terrestris*) (GCF_910591885.1) and the SRR from the WGS of this species (SRR3693403) used in the assigment #8

1. How many reads did not align with the reference genome?

```bash
samtools view -c -f 4 ${BAM}
```
```
107185
```
107,185 reads did not align to the reference genome.

2. How many primary, secondary, and supplementary alignments are in the BAM file?
* primary aligments 
```bash
samtools view -c -F 256 ${BAM}
```
```bash
38156820
```

* secondary aligments 
``` bash 
samtools view -c -f 256 ${BAM}
```
```
0
```

* Suplementary aligments 

``` bash
samtools view -c -f 2048 ${BAM}
```
```
190798
```
using the `flags` `256` and `2048` we can filtrate for the primary, secondary and complementary aligments 


3. How many properly paired alignments on the reverse strand are formed by reads contained in the first pair (read1) file?
```bash
samtools view -c -f 99 ${BAM}
```
```
8681628
```
There are 868,1628 properly paired alignments on the reverse strand in the read 1 file

4. Make a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10.

```bash
samtools view -h -b -q 10 -f 1 -F 256 ${BAM} > ${Q10}
samtools index bam/b_terrestris_q10.bam
```


5. Compare the flagstats for your original and your filtered BAM file.

| **Metric**                                  | **SRA BAM**            | ** Quality 10  BAM**           |
|---------------------------------------------|-----------------------------|----------------------------|
| **Total reads**                             | 38,156,820                  | 28,311,270                 |
| **Primary reads**                           | 37,966,022                  | 28,261,551                 |
| **Secondary reads**                         | 0                           | 0                          |
| **Supplementary reads**                     | 190,798                     | 49,719                     |
| **Duplicates**                              | 0                           | 0                          |
| **Primary duplicates**                      | 0                           | 0                          |
| **Mapped reads**                            | 38,049,635 (99.72%)         | 28,311,270 (100.00%)       |
| **Primary mapped reads**                    | 37,858,837 (99.72%)         | 28,261,551 (100.00%)       |
| **Paired reads**                            | 37,966,022                  | 28,261,551                 |
| **Read1**                                   | 18,983,011                  | 14,135,382                 |
| **Read2**                                   | 18,983,011                  | 14,126,169                 |
| **Properly paired reads**                   | 34,584,800 (91.09%)         | 27,384,193 (96.90%)        |
| **Reads with both mates mapped**            | 37,850,570                  | 28,257,496                 |
| **Singletons (only one mate mapped)**       | 8,267 (0.02%)               | 4,055 (0.01%)              |
| **Reads with mate mapped to different chr** | 762,916                     | 150,766                    |
| **Inter-chromosomal reads (mapQ >= 5)**     | 203,701                     | 150,766                    |

There was a loss of about 25% in the filtered data, but most of the retained reads were primary reads, and low-quality reads were removed.
The percentage of properly paired reads increased from 91.09% to 96.90% in the filtered BAM, suggesting better alignment quality.

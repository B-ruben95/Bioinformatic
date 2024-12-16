# Week 14: Perform a differential expression analysis 
### Ruben Martin

For this assignment, I used the reference genome of the European buff-tailed bumblebee (*Bombus terrestris*) (GCF_910591885.1) and the project [PRJNA355856](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA355856), which uses mRNA-Seq to study the differential gene expression responsible for caste determination in [*Bombus terrestris*]( https://doi.org/10.1111/mec.15752).

## Makefile

The Makefile consists of the following targets:
* `genome`: Downloads the reference genome and the GTF file from NCBI.
* `index`: Indexes the reference genome.
* `design`: Downloads the information from the bioproject and generates the design file.
* `fastq`: Downloads the sequencing reads (SRR).
* `align`: Aligns the RNA reads to the reference genome.
* `count`: Generates the count matrix.
* `csv`: Merges all the counts into a single file.

## Instructions

First, execute the following command to download the genome, index it, and generate the design file:

```bash
make genome index design
```

The Makefile includes a target `rnaseq` that executes both `fastq` and `align`. This target can be run in parallel to automate the download of multiple read files and align them all using the following command:

```bash
cat design.csv | parallel --lb -j 4 --colsep , --header : make rnaseq SRA={run_accession} SAMPLE={sample_alias}
```
To generate the count matrix use the following command:

```bash
make count csv
```


## differentially expressed genes analysis using `DESeq2`

using the target `diff` is posible to make a select differentially expressed genes analysis using `DESeq2`
```bash
make diff
```

```
Input: 13250 rows
# Removed: 1510 rows
# Fitted: 11740 rows
# Significant PVal: 2482 ( 21.1 %)
# Significant FDRs: 1167 ( 9.9 %)
```
2482 genes were found to have a p-value lower than 0.5. After applying a correction factor, 1167 genes were found to be significant.



if this dont work, in the repository there is also a [R script](https://github.com/B-ruben95/Bioinformatic/blob/main/HW14/diff.R) that can be use to make a select differentially expressed genes analysis using `DESeq2` and generate a PCA plot and a heatmap plot

### PCA
![](https://github.com/B-ruben95/Bioinformatic/blob/main/HW14/Image/pca_plot.jpg?raw=true)

The variance explained by the first two principal components is 67% for PC1 and 14% for PC2. This indicates that these two components together capture 81% of the total variability in the data.

The late queen larvae (LQ) and late worker larvae (LW) exhibit similar gene expression patterns, as they cluster closely together in the PCA plot. In contrast, the other larval groups—early queen and worker larvae and mid-stage (median) queen and worker larvae are more dispersed across the plot. This dispersion suggests that these groups have distinct characteristics and more variable gene expression compared to the late-stage larvae.

### heatmap
![](https://github.com/B-ruben95/Bioinformatic/blob/main/HW14/Image/heatmap.jpg?raw=true)

The heat map shows that early queen and worker larvae have similar gene expression patterns, while mid-stage and late queen and worker larvae also share similar gene expression. Additionally, there are differences in the genes expressed during the early and late stages of larval development. This suggests that there is no significant difference in gene expression between the larval development of the different bumblebee castes.


# Week 12: Automating a VCF Calling Pipeline
**Ruben Martin**

For this assignment, I used the reference genome for the European buff-tailed bumblebee (*Bombus terrestris*) (GCF_910591885.1) and the project [PRJNA326162](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA326162/), which utilizes WGS to establish the genetic basis of color polymorphism in this [species](https://www.nature.com/articles/s41598-021-87194-y). 

Here, you can find the [Makefile](https://github.com/B-ruben95/Bioinformatic/blob/main/HW12/Makefile) containing all the commands to call the variants for this project.

## Extract the SRR for the Run

To execute the Makefile, we first need to create the design file with the SRR accession numbers:

```bash
bio search PRJNA326162 -H --csv > design.csv
```
This command downloads the project's information and creates the design file.

## Makefile
The Makefile consists of the following targets:
* `download`: Downloads the reference genome from NCBI.
* `index`: Indexes the reference genome.
* `fastq`: Downloads the runs (SRR).
* `align`: Aligns the runs to the reference genome.
* `variant`: Calls the variants in the runs.
* `merge`: Merges all the variants into a single file.
* `stats`: Generates statistics for the merged file.

The target `all` executes only `download`, `index`, `fastq`, `align`, and `variant` because it is necessary to first generate all the VCF files for the samples before merging them into one file.

## Run the Makefile Using GNU Parallel
The Makefile can be executed in parallel using the following command:

```bash
cat design.csv | head -8 |     parallel --lb -j 4 --colsep , --header :     make all SRA={run_accession} SAMPLE={sample_alias}
```
This project contains 32 SRA files. We will execute the command for a subset of the samples. In this case, the command will process the first 8 samples:

```bash
make all SRR=SRR3692915 SAMPLE=I-D1
make all SRR=SRR3692989 SAMPLE=I-D2
make all SRR=SRR3692994 SAMPLE=I-D3
make all SRR=SRR3692998 SAMPLE=I-D6
make all SRR=SRR3693012 SAMPLE=I-D11
make all SRR=SRR3693014 SAMPLE=I-D12
make all SRR=SRR3693010 SAMPLE=I-D9
```

## Visualize the VCF

To visualize the results, it is first necessary to merge the VCF files into one. Use the following command:

```bash
make merge
```
This will execute the `merge` target to create a single VCF file from all the runs:

```bash
merge:
    # Merge VCF files into a single one.
    bcftools merge -0 vcf/*.vcf.gz -O z > merged.vcf.gz

    # Index the merged VCF file
    bcftools index merged.vcf.gz
```
Once the files are merged, the resulting file can be visualized in `IGV`.

![IGV](https://github.com/B-ruben95/Bioinformatic/blob/main/HW12/Image/Screenshot%202024-11-30%20at%2012.27.34.png?raw=true)

Using the `stats` target, data from the variants can be summarized:

```bash
make stats
```
* **Samples**: 7
* **SNPs**: 1,271,271
* **Indels**: 72,265
* **Total Variants**: 1,343,536
* **Variants with a quality score greater than 30**: 935,093

For this subset of the sequences, there are a total of 1,343,536 variants across 7 samples. Among the variant types, 1,271,271 are SNPs, and 72,265 are indels. This indicates a higher frequency of SNPs compared to indels, which is typical in most genomic datasets. Additionally, 935,093 variants have a quality score greater than 30, while approximately 30% of the variants have low-quality reads.

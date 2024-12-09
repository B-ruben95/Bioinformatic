# Week 13: Generate an RNA-Seq Count Matrix

### Ruben Martin

For this assignment, I used the reference genome of the European buff-tailed bumblebee (*Bombus terrestris*) (GCF_910591885.1) and the project [PRJNA1036000](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1036000/), which uses mRNA-Seq to study the differential gene expression responsible for caste determination in [*Bombus terrestris*](https://link.springer.com/article/10.1007/s13592-024-01117-0).

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
cat design.csv | head -8 | parallel --lb -j 4 --colsep , --header : make rnaseq SRA={run_accession} SAMPLE={sample_alias}
```

You can visualize the `.bw` files in IGV to see the genes expressed in each sample:

![IGV](https://github.com/B-ruben95/Bioinformatic/blob/main/HW13/Image/Screenshot%202024-12-08%20at%2017.40.10.png?raw=true)

Finally, execute the following command to generate the count matrix:

```bash
make count csv
```

Here is a subset of the count matrix:

| Gene          | b_terrestris_B_100 | b_terrestris_B_246 | b_terrestris_B_36 | b_terrestris_B_39 | b_terrestris_B_415 | b_terrestris_B_9 | b_terrestris_G_87 |
|---------------|---------------------|---------------------|--------------------|--------------------|---------------------|------------------|-------------------|
| LOC100648268  | 2                   | 0                   | 0                  | 0                  | 0                   | 0                | 0                 |
| LOC100650023  | 139                 | 138                 | 20                 | 22                 | 102                 | 116              | 6                 |
| LOC125384698  | 0                   | 2                   | 0                  | 0                  | 0                   | 0                | 0                 |
| LOC100650139  | 142                 | 115                 | 52                 | 56                 | 91                  | 115              | 66                |
| LOC125384718  | 1                   | 0                   | 0                  | 0                  | 1                   | 6                | 1                 |
| LOC100650265  | 92                  | 96                  | 85                 | 65                 | 73                  | 92               | 71                |
| LOC100650388  | 24                  | 22                  | 28                 | 14                 | 29                  | 91               | 26                |
| LOC100650511  | 77                  | 78                  | 21                 | 18                 | 112                 | 108              | 124               |
| LOC105666912  | 61                  | 26                  | 44                 | 29                 | 34                  | 41               | 39                |
| LOC100650629  | 52                  | 104                 | 95                 | 44                 | 109                 | 65               | 56                |
| LOC100650754  | 299                 | 277                 | 84                 | 90                 | 284                 | 118              | 244               |
| LOC125384641  | 42                  | 46                  | 31                 | 24                 | 125                 | 128              | 73                |
| LOC100650872  | 205                 | 328                 | 153                | 146                | 300                 | 435              | 214               |
| LOC100651113  | 444                 | 653                 | 429                | 176                | 698                 | 548              | 433               |
| LOC100648385  | 89                  | 68                  | 160                | 203                | 114                 | 302              | 95                |
| LOC125384793  | 0                   | 0                   | 0                  | 0                  | 0                   | 0                | 0                 |
| LOC125385192  | 0                   | 8                   | 0                  | 0                  | 4                   | 0                | 0                 |

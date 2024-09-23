# Week 3: Download, Generate and Visualize GFF files in IGV
**Ruben Martin**

I use the NCBI data base to download the references genome from Escheriquia coli, its one of the model organis in molecular biologi and its genome is pretty well anotated

To view the information of the file before downlading i use the folowing command line
```bash
datasets summary genome accession GCF_000005845.2 | jq
```
and it print the following:
```
New version of client (16.28.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets.
{
  "reports": [
    {
      "accession": "GCF_000005845.2",
      "annotation_info": {
        "name": "Annotation submitted by NCBI RefSeq",
        "provider": "NCBI RefSeq",
        "release_date": "2022-03-09",
        "stats": {
          "gene_counts": {
            "non_coding": 206,
            "protein_coding": 4288,
            "pseudogene": 145,
            "total": 4639
          }
        }
      },
      "assembly_info": {
        "assembly_level": "Complete Genome",
        "assembly_name": "ASM584v2",
        "assembly_status": "current",
        "assembly_type": "haploid",
        "bioproject_accession": "PRJNA225",
        "bioproject_lineage": [
          {
            "bioprojects": [
              {
                "accession": "PRJNA225",
                "title": "Model organism for genetics, physiology, biochemistry"
              }
            ]
          }
        ]
  ...
```
with this we can see relevant infortmation about the genome, like accession number, publication date, number if genes, number of coding sequens (CDS), etc. 

to download the genemo i use the following command 
```bash
datasets download genome accession GCF_000005845.2 --include gff3,genome
```
it print

```bash
New version of client (16.28.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    1.82MB valid zip structure -- files not checked
Validating package [================================================] 100% 6/6
```
here we include the genome and the gff3 file

when we upload the files to the 
![screenshot](https://github.com/B-ruben95/Bioinformatic/blob/main/lecture4/IGV_all_infomation.png?raw=true)

we can see that the E coli genome have only one chromosome and has a length of 4641 kb. 

with this line of command we cut only the lines that have gene information to make it more simple to visualize 
```bash
cat E_coli/genomic.gff | awk '$3 == "gene"' > genes.gff
```
![screenshot](https://github.com/B-ruben95/Bioinformatic/blob/main/lecture4/IGV_only_genes.png?raw=true)
here we can visualice only the genes

now i'm going to create my gff3 file using the next file 
```
NC_000913.3	.	gene 	523260	527541	.	+	.	Parent=transcript1;ID=gene1
NC_000913.3	.	gene 	527581	527949	.	+	.	Parent=transcript1;ID=gene2
```
![screensshsot](https://github.com/B-ruben95/Bioinformatic/blob/main/lecture4/IGV_my_gff3.png?raw=true)

when the file have the gene feature it dont generate the fish bone conection, when we changes the feature to CDS it show the fish bone conection showing that the two sequences are related 
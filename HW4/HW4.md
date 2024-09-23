# Week 4. Scripting and Sequence Ontologies
**Ruben Martin**

### part 1 Write a script.
1. Run your script on your original data and verify that it works

```bash
bash hw4.sh
```

```bash
+ NCBI=GCF_000005845.2
+ TYPE=gff3,genome
+ GENOME='ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2*.fna'
+ GENOMIC=ncbi_dataset/data/GCF_000005845.2/genomic.gff
+ '[' -f ncbi_dataset.zip ']'
+ datasets download genome accession GCF_000005845.2 --include gff3,genome
New version of client (16.29.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    1.82MB valid zip structure -- files not checked
Validating package [================================================] 100% 6/6
+ '[' -d ncbi_dataset ']'
+ unzip ncbi_dataset.zip
Archive:  ncbi_dataset.zip
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna  
  inflating: ncbi_dataset/data/GCF_000005845.2/genomic.gff  
  inflating: ncbi_dataset/data/dataset_catalog.json  
  inflating: md5sum.txt              
+ samtools faidx ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna
+ cat 'ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2*.fna.fai' ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna.fai
+ head
NC_000913.3	4641652	72	80	81
+ seqkit stats ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna
file                                                                    format  type  num_seqs    sum_len    min_len    avg_len    max_len
ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna  FASTA   DNA          1  4,641,652  4,641,652  4,641,652  4,641,652
+ cat ncbi_dataset/data/GCF_000005845.2/genomic.gff
+ awk '$3 == "gene"'
+ echo 'the GFF3 file has been created with only the genes'
the GFF3 file has been created with only the genes
+ cat genes.gff
+ wc -l
    4494
```
my script run successfully 

2. 

### Part 2: Make use of ontologies


1. Choose a feature type from the GFF file and look up its definition in the sequence ontology.

```bash
bio explain cds
```

print 

```bash
## cds (SO:0000316)

A contiguous sequence which begins with, and includes, a
start codon and ends with, and includes, a stop codon.
```

2. Find both the parent terms and children nodes of the term

```
Parents:
- mrna_region 

Children:
- polypeptide (derives_from)
- cds_region (part_of)
- edited_cds 
- cds_fragment 
- transposable_element_cds 
- cds_extension 
- cds_independently_known 
- cds_predicted 
```

3. Provide a short discussion of what you found.
CDS 




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
+ cat genes.gff
+ wc -l
    4494
```

2. run your script on their data. If the script is reusable, you can replace your variables with theirs and run the script

Noah Alexander Kinscherf use the ensembl web page to downlaod the genome, and in my scipt i use the NCBI, so i'll have to modify the variables and the commands to by able to run their data in my script. in this case my script is not reusable becauses of the structure of the commands. 

my script use the number of accession to downlaod the genome and GFF3 file.

```bash
#variables

NCBI=GCF_013318015.2
TYPE=gff3,genome
GENOME=ncbi_dataset/data/${NCBI}/${NCBI}*.fna
GENOMIC=ncbi_dataset/data/$NCBI/genomic.gff
```
and noah use a link to the ensembl repository 
```bash
GENOME="Ursus_maritimus.UrsMar_1.0.112.gff3"
URL=https://ftp.ensembl.org/pub/current_gff3/ursus_maritimus/Ursus_maritimus.UrsMar_1.0.112.gff3.gz
```

in the caises of Nikol Chantzi she did use the NCBI to download the genemo of the Human in this cases my script is reusable, i replace de accession number of E. coli  with the one from Human in my script

```bash
NCBI=GCF_013318015.2
TYPE=gff3,genome
GENOME=ncbi_dataset/data/${NCBI}/${NCBI}*.fna
GENOMIC=ncbi_dataset/data/$NCBI/genomic.gff
```
and run the scrip

```bash
bash hw4_val.sh
```
```bash
NCBI=GCF_000001405.40
+ TYPE=gff3,genome
+ GENOME='ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40*.fna'
+ GENOMIC=ncbi_dataset/data/GCF_000001405.40/genomic.gff
+ '[' -f ncbi_dataset.zip ']'
+ datasets download genome accession GCF_000001405.40 --include gff3,genome
New version of client (16.29.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    1.05GB valid zip structure -- files not checked
Validating package [================================================] 100% 6/6
+ '[' -d ncbi_dataset ']'
+ unzip ncbi_dataset.zip
Archive:  ncbi_dataset.zip
  inflating: README.md               
  inflating: ncbi_dataset/data/assembly_data_report.jsonl  
  inflating: ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna  

  inflating: ncbi_dataset/data/GCF_000001405.40/genomic.gff  
  inflating: ncbi_dataset/data/dataset_catalog.json  
  inflating: md5sum.txt              
+ samtools faidx ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
+ cat 'ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40*.fna.fai' ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna.fai
+ head
NC_000001.11    248956422       69      80      81
NT_187361.1     175055  252068568       80      81
NT_187362.1     32032   252245933       80      81
NT_187363.1     127682  252278487       80      81
NT_187364.1     66860   252407887       80      81
NT_187365.1     40176   252475704       80      81
NT_187366.1     42210   252516504       80      81
NT_187367.1     176043  252559363       80      81
NT_187368.1     40745   252737728       80      81
NT_187369.1     41717   252779104       80      81
+ seqkit stats ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna
file                                                                        format  type  num_seqs        sum_len  min_len    avg_len      max_len
ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna  FASTA   DNA        705  3,298,430,636      970  4,678,625  248,956,422
+ cat ncbi_dataset/data/GCF_000001405.40/genomic.gff
+ awk '$3 == "gene"'
+ echo 'the GFF3 file has been created with only the genes'
the GFF3 file has been created with only the genes
+ cat genes.gff
+ wc -l
   47532
   ```

For this data set, my script is reusable

3. Add more functions to the script that also print some of their results. Were you able to reproduce their results? Make a note in the report

* nicole 
```
bash feature_counter.sh GCF_000001405.40
```
when you run Nikole script this what it print if you provite and accession number

```bash
Detected accession GCF_000001405.40.
Processing accession GCF_000001405.40.
Red Alert! Red Alert!
Unable to detect accession. Initializing download protocol 4032

New version of client (16.29.0) available at https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/mac/datasets.
Collecting 1 genome record [================================================] 100% 1/1
Downloading: ncbi_dataset.zip    77.8MB valid zip structure -- files not checked
Validating package [================================================] 100% 5/5
GFF File has been succesfully downloaded for accession GCF_000001405.40.
Food is in the oven. Please wait...
Total features detected: 4903918
Process has completed. Your resuls are ready here: GCF_000001405.40.counts.txt!
```
i was able to recreate her results 

* noah
```bash
bash week4.sh 
```
print

```bash
The number of sequence regions is    23819

435225	exon
415739	CDS
142264	biological_region
34027	mRNA
23819	region
18724	gene
15455	five_prime_UTR
11224	three_prime_UTR
5296	ncRNA_gene
3384	lnc_RNA

The number of features of type gene in this assembly is    18724

The number of features in the document is  1108663
```


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







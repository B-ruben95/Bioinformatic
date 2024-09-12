1. Tell us a bit about the organism.
The file I chose was "Ailuropoda_melanoleuca.ASM200744v2.112.gff3.gz," which corresponds to the giant panda. This mammal is endemic to China and is one of the most endangered species of mammals in the world due to its limited geographical range and deforestation.

2. How many features does the file contain?

```bash
cat Ailuropoda_melanoleuca.ASM200744v2.112.gff3 | grep -v "^#" | wc -l
```

 print
```bash
1283346
```
The gff3 file of the the panda have 128 3346 features. 

3.How many sequence regions (chromosomes) does the file contain?
    cat Ailuropoda_melanoleuca.ASM200744v2.112.gff3 | grep sequence-region | wc -l
    77287

the panda has a highly fragmented genome and does not have its genome organized by chromosomes. According to the file, it contains 77,287 sequence regions.

4. How many genes are listed for this organism?

    cat panda.gff3 | cut -f 3| sort | uniq -c | sort -rn | head                  
    451506 exon
    430935 CDS
    202966 biological_region
    77287 region
    39273 mRNA
    31558 five_prime_UTR
    21284 three_prime_UTR
    20857 gene
    3022 ncRNA_gene
    1289 snRNA

This organism have 20 857 genes

5. What are the top-ten most annotated feature types (column 3) across the genome?

    cat panda.gff3 | cut -f 3| sort | uniq -c | sort -rn | head                  
    451506 exon
    430935 CDS
    202966 biological_region
    77287 region
    39273 mRNA
    31558 five_prime_UTR
    21284 three_prime_UTR
    20857 gene
    3022 ncRNA_gene
    1289 snRNA

The most annotated parts of the genome are exons, CDS, biological regions, regions, mRNA, 5' UTR, 3' UTR, genes, ncRNA's and scRNA's

6.Having analyzed this GFF file, does it seem like a complete and well-annotated organism?

Despite being one of the first mammals to have its genome sequenced, the organization of the genome is quite poorly done. It lacks chromosome and remains highly fragmented. Nevertheless, a significant number of genes, RNA, exons, and other elements are mapped within the genome.

'''bash

'''
# Week 10: Generating a Variant Call File
**Ruben Martin**

The [makefile](https://) contains all the commands to call variants from the SRA dataset SRR3693403 of *Bombus terrestris*.

## Variant Calling

```bash
gunzip -c vcf/b_terrestris.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq -c
```
```bash
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  bam/b_terrestris.bam
NC_063269.1     31388   .       A       C       4.38466 .       DP=1;AD=0,1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=32    GT:PL:DP:SP:ADF:ADR:AD:GQ       1/1:32,3,0:1:0:0,1:0,0:0,1:127
NC_063269.1     31403   .       T       C       4.38466 .       DP=1;AD=0,1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,1,0;MQ=32    GT:PL:DP:SP:ADF:ADR:AD:GQ       1/1:32,3,0:1:0:0,1:0,0:0,1:127
NC_063269.1     40691   .       G       T       37.2985 .       DP=14;AD=11,3;VDB=0.460446;SGB=-0.511536;RPBZ=-1.47915;MQBZ=3.22351;MQSBZ=-3.54897;BQBZ=-0.894167;SCBZ=0;MQ0F=0.714286;AC=2;AN=2;DP4=1,10,3,0;MQ=8      GT:PL:DP:SP:ADF:ADR:AD:GQ        1/1:65,9,0:14:20:1,3:10,0:11,3:127
NC_063269.1     40847   .       G       A       4.38466 .       DP=1;AD=0,1;SGB=-0.379885;MQ0F=0;AC=2;AN=2;DP4=0,0,0,1;MQ=32    GT:PL:DP:SP:ADF:ADR:AD:GQ       1/1:32,3,0:1:0:0,0:0,1:0,1:127
```
To call the variants, I used the following command:

```bash
echo "SNPs: $(grep "SN[[:space:]]0[[:space:]]number of SNPs" stats_report.txt | awk '{print $6}')"
echo "Indels: $(grep "SN[[:space:]]0[[:space:]]number of indels" stats_report.txt | awk '{print $6}')"
echo "Total Variants: $(grep "SN[[:space:]]0[[:space:]]number of records" stats_report.txt | awk '{print $6}')"
```

### Results
```bash
Total Variants: 1130592
SNPs: 1054056
Indels: 76536
Variants with a quality score greater than 30: 871577
```

In total, there are `1,130,592` variants in this SRA dataset from the *Bombus terrestris* genome. Among these, I identified `1,054,256` SNPs and `76,536` indels. This indicates a higher frequency of SNPs compared to indels, which is common in most genomic datasets. Additionally, `871,577` variants have a quality score greater than 30, meaning approximately 30% of the variants have low-quality reads.

### Variant Counts Per Chromosome
```bash
# Variants per chromosome
gunzip -c ${VCF} | grep -v "^#" | cut -f1 | sort | uniq -c
# SNPs per chromosome
gunzip -c ${VCF} | grep -v "^#" | awk 'length($4) == 1 && length($5) == 1' | cut -f1 | sort | uniq -c
# Indels per chromosome
gunzip -c ${VCF} | grep -v "^#" | awk 'length($4) != 1 || length($5) != 1' | cut -f1 | sort | uniq -c
```

I also generated a table showing the number of variants by chromosome. For clarity, only Ensembl chromosomes are displayed:

| Chromosome   | Count  |  SNP  | Indel |
|--------------|--------|-------|-------|
| NC_063269.1  | 66398  | 61670 | 4728  |
| NC_063270.1  | 65624  | 60684 | 4940  |
| NC_063271.1  | 77178  | 71507 | 5671  |
| NC_063272.1  | 62605  | 57664 | 4941  |
| NC_063273.1  | 46051  | 42528 | 3523  |
| NC_063274.1  | 77804  | 72037 | 5767  |
| NC_063275.1  | 90234  | 83875 | 6359  |
| NC_063276.1  | 45620  | 42085 | 3535  |
| NC_063277.1  | 70757  | 65656 | 5101  |
| NC_063278.1  | 74864  | 69426 | 5438  |
| NC_063279.1  | 76826  | 71110 | 5716  |
| NC_063280.1  | 49118  | 45256 | 3862  |
| NC_063281.1  | 55597  | 51612 | 3985  |
| NC_063282.1  | 42669  | 39355 | 3314  |
| NC_063283.1  | 43070  | 39665 | 3405  |
| NC_063284.1  | 34438  | 31914 | 2524  |
| NC_063285.1  | 38491  | 36119 | 2372  |
| NC_063286.1  | 15284  | 14277 | 1007  |

![IGV Screenshot](https://github.com/B-ruben95/Bioinformatic/blob/main/HW10/image/Screenshot%202024-11-10%20at%2016.14.21.png?raw=true)

## Identifying Issues with the Variant Caller: False Positives and False Negatives
![IGV Screenshot](https://github.com/B-ruben95/Bioinformatic/blob/main/HW10/image/Screenshot%202024-11-10%20at%2022.53.42.png?raw=true)
![IGV Screenshot](https://github.com/B-ruben95/Bioinformatic/blob/main/HW10/image/Screenshot%202024-11-10%20at%2022.53.51.png?raw=true)

The ends of the chromosomes have poor read coverage, which may lead to incorrect variant calls, such as false positives and false negatives.

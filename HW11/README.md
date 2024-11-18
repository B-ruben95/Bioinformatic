
# Week 11: Establish the Effects of Variants

**Ruben Martin**

For this assignment, I used the reference genome of the European buff-tailed bumblebee (Bombus terrestris) (GCF_910591885.1) and the SRR from the WGS of this species (SRR3693403).

Here is the [Makefile](https://github.com/B-ruben95/Bioinformatic/blob/main/HW11/Makefile#:~:text=Image-,Makefile,-HW2) containing all the commands to establish the effects of variants.

| Statistic                        | Value     |
|-----------------------------------|-----------|
| Lines of input read               | 1,130,592 |
| Variants processed                | 1,130,592 |
| Variants filtered out             | 0         |
| Novel / existing variants         | -         |
| Overlapped genes                  | 13,124    |
| Overlapped transcripts            | 31,034    |
| Overlapped regulatory features    | -         |

VEP generated the following graphs:

### Variant Classes

![Variant classes](https://github.com/B-ruben95/Bioinformatic/blob/main/HW11/Image/Screenshot%202024-11-17%20at%2016.14.14.png?raw=true)

| Variant Class         | Count     |
|-----------------------|-----------|
| sequence_alteration   | 574       |
| insertion             | 36,898    |
| deletion              | 39,064    |
| SNV                   | 1,054,056 |

The analysis identified a total of 1,054,056 single nucleotide variants (SNVs), along with 36,898 insertions and 39,064 deletions.

### Most Severe Consequences

![Consequences (most severe)](https://github.com/B-ruben95/Bioinformatic/blob/main/HW11/Image/Screenshot%202024-11-17%20at%2016.14.06.png?raw=true)

| Consequence Type                      | Count    |
|---------------------------------------|----------|
| intron_variant                        | 665,423  |
| intergenic_variant                    | 268,159  |
| 3_prime_UTR_variant                   | 35,444   |
| non_coding_transcript_exon_variant    | 23,254   |
| 5_prime_UTR_variant                   | 15,908   |
| synonymous_variant                    | 15,116   |
| missense_variant                      | 14,736   |
| upstream_gene_variant                 | 56,241   |
| downstream_gene_variant               | 29,952   |
| splice_region_variant                 | 2,118    |
| splice_polypyrimidine_tract_variant   | 1,895    |
| splice_donor_region_variant           | 401      |
| splice_donor_5th_base_variant         | 166      |
| inframe_deletion                      | 108      |
| inframe_insertion                     | 104      |
| stop_gained                           | 945      |
| splice_donor_variant                  | 241      |
| splice_acceptor_variant               | 216      |
| frameshift_variant                    | 64       |
| start_lost                            | 42       |
| stop_lost                             | 39       |
| stop_retained_variant                 | 16       |
| protein_altering_variant              | 2        |
| coding_sequence_variant               | 2        |

Among the most severe variants, intron variants (665,423) were the most abundant. These mutations are usually harmless, but if an intron contains regulatory sequences like splicing sites or elements that affect protein translation, it can have significant effects.

### All Consequences

![Consequences (all)](https://github.com/B-ruben95/Bioinformatic/blob/main/HW11/Image/Screenshot%202024-11-17%20at%2016.14.22.png?raw=true)

| Consequence Type                      | Count       |
|---------------------------------------|-------------|
| intron_variant                        | 3,124,187   |
| downstream_gene_variant               | 481,770     |
| upstream_gene_variant                 | 468,937     |
| non_coding_transcript_variant         | 317,734     |
| intergenic_variant                    | 268,732     |
| 3_prime_UTR_variant                   | 124,608     |
| 5_prime_UTR_variant                   | 32,641      |
| non_coding_transcript_exon_variant    | 52,552      |
| synonymous_variant                    | 50,999      |
| missense_variant                      | 40,161      |
| splice_polypyrimidine_tract_variant   | 10,192      |
| splice_region_variant                 | 7,826       |
| stop_gained                           | 2,460       |
| splice_donor_region_variant           | 1,199       |
| splice_donor_variant                  | 616         |
| splice_acceptor_variant               | 577         |
| splice_donor_5th_base_variant         | 482         |
| inframe_insertion                     | 379         |
| inframe_deletion                      | 360         |
| frameshift_variant                    | 123         |
| start_lost                            | 63          |
| stop_lost                             | 65          |
| stop_retained_variant                 | 31          |
| start_retained_variant                | 9           |
| coding_sequence_variant               | 7           |
| protein_altering_variant              | 5           |

### Coding Consequences

![Coding consequences](https://github.com/B-ruben95/Bioinformatic/blob/main/HW11/Image/Screenshot%202024-11-17%20at%2016.14.30.png?raw=true)

| Consequence Type        | Count    |
|-------------------------|----------|
| synonymous_variant      | 50,999   |
| missense_variant        | 40,161   |
| stop_gained             | 2,460    |
| inframe_insertion       | 379      |
| inframe_deletion        | 360      |
| frameshift_variant      | 123      |
| stop_lost               | 65       |
| start_lost              | 63       |
| stop_retained_variant   | 31       |
| start_retained_variant  | 9        |
| coding_sequence_variant | 7        |
| protein_altering_variant| 5        |

Most variants in the coding regions are synonymous, meaning they do not alter the translated protein. The next most common variants are missense mutations, which can potentially change protein function. Nonsense mutations, the third most abundant, introduce stop codons that can produce truncated, often nonfunctional and deleterious proteins.

### Variants by Chromosome

| Chromosome   | Count  |
|--------------|--------|
| NC_063269.1  | 66,398 |
| NC_063270.1  | 65,624 |
| NC_063271.1  | 77,178 |
| NC_063272.1  | 62,605 |
| NC_063273.1  | 46,051 |
| NC_063274.1  | 77,804 |
| NC_063275.1  | 90,234 |
| NC_063276.1  | 45,620 |
| NC_063277.1  | 70,757 |
| NC_063278.1  | 74,864 |
| NC_063279.1  | 76,826 |
| NC_063280.1  | 49,118 |
| NC_063281.1  | 55,597 |
| NC_063282.1  | 42,669 |
| NC_063283.1  | 43,070 |
| NC_063284.1  | 34,438 |
| NC_063285.1  | 38,491 |
| NC_063286.1  | 15,284 |
| scaffold     | 289    |

Some chromosomes exhibited a higher density of variants.

# Week 6: FASTQ Quality Control
**Ruben Martin**

I selected the whole genome data for Crabro peltarius with the accession number SRR12432506. This is a wasp native to Europe.
The sequences are part of the BioProject PRJNA607895 focusing on Danish organisms. This project contains 1,165 SRA experiments of different taxonomic groups.
A paper uses data from this project, but it only includes vertebrates in its analysis (doi: 10.1002/edn3.138). Currently, there are no publications that evaluate the data for other taxa like insects.

1. Evaluation of the Quality of the Downloaded Data
The sequencing quality is not very high, as the reads contain Illumina adapters. Below are the quality control results for Read 1 and Read 2.

Read 1 Quality

![Q](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Q.jpg?raw=true)

Adapters in the Data
![adapters](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Quality_Read_Adapters.jpg?raw=true)


Read 2 Quality

![Q2](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Quality_Read2.png?raw=true)


Adapters Detected
![adapter](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Screenshot%202024-10-07%20at%2000.40.27.png?raw=true)

2. Improving the Quality of the Reads in the Dataset
After quality improvement, Read 1 shows better results than Read 2. However, Read 2 requires additional filtering to achieve higher quality, though this would lead to some data loss, especially in Read 1.

Both reads still require better trimming at the beginning of the sequences. Bases 1 to 6 contain some errors.

Trimmed Read 1

![adapter](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Trimmed1.png?raw=true)


![adapter](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Screenshot%202024-10-07%20at%2000.43.51.png?raw=true)


Trimmed Read 2

![adapter](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/trimmed2.png?raw=true)



![adapter](https://github.com/B-ruben95/Bioinformatic/blob/main/HW6/Image/Screenshot%202024-10-07%20at%2000.43.51.png?raw=true)


Acording to the NCBI page for this project, this dataset was sequenced for the purpose of assembling genomes for future reference. The authors need to perform additional cleaning on the files to improve data quality.

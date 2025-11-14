# Motif Base Editing Screen

### Overview:
This repository includes the pipeline for analysing NGS count data from the Motif Dependency Map described in: 

**“A proteome-wide dependency map of protein interaction motifs”** Ambjørn SM, Meeusen B,Kliche J, Wang L, Garvanska DH, Kruse T, Mendez BL, Mann M, Mailand N, Hertz EPT*, Davey NE*, Nilsson J*. 
*bioRxiv DOI:2024.09.11.612445;*

### Pipeline for processing motif screen base editing data

#### Prerequisites:


The code base is containerised in a Docker container and therefore requires that docker is installed on your system. Information on installing Docker is available at the Docker website [here](https://docs.docker.com/get-started/get-docker/).


The pipelines requires the EnsEMBL Homo_sapiens GRCh38 in fasta format. The file is available here: 

```
https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Unzip the fasta file in the ./be_annotation_data/genomes/ folder. 


#### Running the pipeline:

Once Docker is running the container can be built, run and accessed. First, build and run the docker image using the following command:
```
docker-compose -f be_docker.yml up --detach 
# or depending on version
docker compose -f be_docker.yml up --detach 
```

Once the docker file is running the docker file can be accessed using the following command (Docker will download the required software on the first run which may take up to 10 minutes):
```
docker exec -it be_docker-1 bash
# if be_docker-1 container does not exists use the command 'docker container ls' to find the name of the container on your system
```
Docker will start an Ubuntu container with the required dependencies preinstalled with an open bash. 

If running for the first time, the human genome fasta file must be processed with makeblastdb before running the pipeline.
```
makeblastdb -dbtype nucl  -in "/home/be_annotation_data/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa" -out GRCh38_dna_primary_assembly -parse_seqids
```

The processing pipeline can then be run using supplied runBaseEditingScreenParser.py
```
#Test by running
python run_be_screen_pipeline.py --job_file ./job_files/motif/P53.guides.test.json
# This will pull data for one protein and should complete if all parts of the pipeline are correctly installed.

#Test by running 
python run_be_screen_pipeline.py --job_file ./job_files/motif/ABE_NGN.json 
```

There are four jobfile, one for each subscreen:
```
./job_files/motif/ABE_NGN.json
./job_files/motif/ABE_NGG.json
./job_files/motif/CBE_NGN.json
./job_files/motif/CBE_NGG.json
```

This pipeline will perform several tasks including processing the NGS count data, annotating the gRNAs and motif regions, and calculating enrichment statistics. The scripts download a large amount of data and require an active internet connection. 

The data fo gRNA is processed to define significantly changing gRNAs using a limma statistical analysis . Next generation sequencing data from the day 0 and 18 samples are converted to gRNA sequencing counts by mapping to motif-targeting and control gRNAs (essential splice sites, non-targeting , intergenic) in the motif dependency map gRNA library design. Three replicates were collected per time point, resulting in three data points per gRNA per time point. 

Low-abundance gRNAs are removed (counts <30). The gRNAs are mapped to the genome to find non-specific gRNAs with off-target genome binding sites. Then gRNAs with 5 or more off-target matches in the genome were removed from downstream analysis. 

Raw sequencing counts are normalised per sample to log2 transcripts per million (log2TPM). The log2TPM values were compared using limma, resulting in fold change (FC) and p-values for each gRNA. These p-values were compared to the positive and negative controls, and a p-value cut-off of 0.01 showed a false positive rate of ~1%. Consequently, a p-value of less than 0.01 was set as the threshold considered to have an effect on cell proliferation. 

Next, the gRNAs mapping to each motif-containing region were analysed to define groups of significantly changing gRNAs. The main metrics used was the number of unique significant gRNAs per motif across all screens. Per screen, the direction of the effect on proliferation was quantified as the mean fold change of the gRNAs to define depletion or enrichment, and the fold change of the gRNAs for a region were compared against the fold change for all gRNAs in the analysis using a Mann–Whitney U test to calculate an enrichment or depletion p-value. 

After ~10 minutes the results will appear in 
```
./be_results/
```

*Tested with:*

Hardware/Operating system: Apple M2 Pro, 32 GB RAM, Docker version 24.0.6 running Ubuntu 22.04

Dependency versions: biopython 1.81, seaborn 0.11.2, matplotlib 3.7.5, numpy 1.24.4, scipy 1.10.1, pandas 2.0.3, rpy2 3.5.4 in python 3.7.17 ; limma 3.60.6 in R version 3.6.3; 
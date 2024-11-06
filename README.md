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

Once the docker file is running the docker file can be accessed using the following command:
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
python run_be_screen_pipeline.py --job_file ./job_files/motif/ABE_NGN.json 
```

There are four jobfile, one for each subscreen:
```
./job_files/motif/ABE_NGN.json
./job_files/motif/ABE_NGG.json
./job_files/motif/CBE_NGN.json
./job_files/motif/CBE_NGG.json
```

The scripts dowload a large amount of data and require an active internet connection.
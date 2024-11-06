# Motif Base Editing Screen
### Tools for processing motif screen base editing data

“A proteome-wide dependency map of protein interaction motifs” Ambjørn SM, Meeusen B,Kliche J, Wang L, Garvanska DH, Kruse T, Mendez BL, Mann M, Mailand N, Hertz EPT*, Davey NE*, Nilsson J*. 
bioRxiv DOI:2024.09.11.612445;

#### Running the pipeline:

The pipelines requires the EnsEMBL Homo_sapiens GRCh38 in fasta format. The is available here: https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz. Unzip the fasta file in the /home/be_annotation_data/genomes/ folder. The fasta file must then processed with makeblastdb before running the pipeline.
```
/home/be_annotation_data/genomes/
makeblastdb -dbtype nucl  -in "/home/be_annotation_data/genomes/Homo_sapiens.GRCh38.dna.primary_assembly.fa" -out GRCh38_dna_primary_assembly -parse_seqids
```

The code base is containerised in a Docker container and therefore requires that docker is installed on your system. Information on installing Docker is avaialble at the Docker website [here](https://docs.docker.com/get-started/get-docker/). Once Docker is running the container can be installed and run using the supplied shell script. This script will build an ubuntu container containing all the dependencies required for the analysis pipeline. 

The Docker is built and run using the following command.

```
sh ./run_base_editing_processing_docker.sh 
```

or directly call the following docker commands

```
docker-compose -f be_docker.yml up --detach 
# or depending on version
docker compose -f be_docker.yml up --detach 
docker exec -it be_docker-1 bash
# if be_docker-1 container does not exists use the command 'docker container ls' to find the name of the container on your system
```

The Docker will start with an open bash. The processing pipeline can then be run using supplied runBaseEditingScreenParser.py
```
python run_be_screen_pipeline.py --job_file ./job_files/motif/ABE_NGN.json 
```

The scripts dowload a large amount of data and require an active internet connection.

There are four jobfile, one for each subscreen:
```
./job_files/motif/ABE_NGN.json
./job_files/motif/ABE_NGG.json
./job_files/motif/CBE_NGN.json
./job_files/motif/CBE_NGG.json
```

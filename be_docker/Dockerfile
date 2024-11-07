FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get update --fix-missing
RUN apt-get install -y wget
RUN apt-get install  -y apt-utils
RUN apt-get -y upgrade

RUN apt-get update --fix-missing
RUN apt-get install -y supervisor
RUN apt-get install -y host net-tools adduser htop lshw
RUN apt-get update --fix-missing
RUN useradd -ms /bin/bash submitter

ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

RUN chmod -R 777 /var/log/supervisor

RUN apt-get update --fix-missing
RUN apt-get install -y python2
RUN apt-get update --fix-missing
RUN apt-get install -y python3-pip python3-dev
RUN ln -s /usr/bin/python3 /usr/local/bin/python
RUN pip3 install --upgrade pip
RUN apt-get update --fix-missing

RUN pip install --upgrade pip
RUN pip install --upgrade setuptools

RUN /usr/bin/python3 -m pip install --upgrade pip
RUN pip install requests

RUN pip install biopython>=1.81
RUN pip install pytest
RUN pip install dataset

RUN pip install futures
RUN pip install numpy
RUN pip install scipy>=1.9.0
RUN pip install scikit-learn
RUN pip install unidecode

RUN pip install matplotlib>=3.5.3
RUN pip install seaborn>=0.11.2
RUN pip install pandas>=1.4.3

RUN pip install networkx
RUN python3 -m pip install --user svg.py
RUN python3 -m pip install --user cairosvg

#### DSSP
#WORKDIR     /usr/bin
#RUN wget https://github.com/PDB-REDO/dssp/releases/download/v4.4.0/mkdssp-4.4.0-linux-x64
#RUN chmod a+x mkdssp-4.4.0-linux-x64
#RUN ln -s mkdssp-4.4.0-linux-x64 mkdssp
RUN apt-get install -y dssp

#### BLAST
RUN apt-get install -y ncbi-blast+

### R
RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base

RUN apt-get install python3.10-dev build-essential gcc
RUN apt-get install -y libpcre2-dev libbz2-dev zlib1g-dev liblzma-dev libicu-dev libblas-dev
RUN pip install --upgrade pip
RUN pip install rpy2>=3.5.13

### Move to correct directory
WORKDIR /home/scripts/

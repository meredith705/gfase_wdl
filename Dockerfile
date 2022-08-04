FROM ubuntu:20.04
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu

RUN mkdir -p /home/apps

RUN cd /home/apps && \
    apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget make build-essential cmake \
    protobuf-compiler pkg-config libprotobuf-dev libjansson-dev libhts-dev libncurses-dev libbz2-dev liblzma-dev


RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && tar xvjf samtools-1.15.1.tar.bz2 && \
    rm samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure --prefix /usr/local/bin && \
    make


RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && tar -xf bwa-mem2-2.2.1_x64-linux.tar.bz2 \
    && mv bwa-mem2-2.2.1_x64-linux /usr/local/bin/


#ADD https://api.github.com/repos/rlorigro/GFAse/git/refs/heads/wdl_debug version.json

RUN cd /home/apps && \
    git clone https://github.com/rlorigro/GFAse.git && \
    cd GFAse && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make -j 8 
    
ENV PATH="/home/apps/GFAse/build:${PATH}"
ENV PATH="/home/apps/GFAse/data:${PATH}"

RUN mkdir -p /data
WORKDIR /data


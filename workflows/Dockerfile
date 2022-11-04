FROM ubuntu:20.04
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu

RUN mkdir -p /home/apps


RUN cd /home/apps && \
    apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget make build-essential cmake \
    python3.8 python3.8-dev python3-pip \
    protobuf-compiler pkg-config libprotobuf-dev libjansson-dev libhts-dev libncurses-dev \
    libbz2-dev liblzma-dev zlib1g-dev autoconf libcurl4-openssl-dev curl libomp-dev

WORKDIR /home/apps
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    tar xvf samtools-1.15.1.tar.bz2 && \
    rm samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure && \
    make

WORKDIR /home/apps
RUN wget https://github.com/bwa-mem2/bwa-mem2/releases/download/v2.2.1/bwa-mem2-2.2.1_x64-linux.tar.bz2 && \
    tar -jxf bwa-mem2-2.2.1_x64-linux.tar.bz2

WORKDIR /home/apps
RUN mkdir kmc && \
    cd kmc && \
    wget https://github.com/refresh-bio/KMC/releases/download/v3.2.1/KMC3.2.1.linux.tar.gz && \
    tar -xzf KMC3.2.1.linux.tar.gz

RUN cd /home/apps && \
    git clone https://github.com/rlorigro/GFAse.git && \
    cd GFAse && \
    git checkout 5412b86db184bed984129ac17d811bf554308afb && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

ENV PATH="/home/apps/GFAse/build:${PATH}"
ENV PATH="/home/apps/GFAse/data:${PATH}"
ENV PATH="/home/apps/GFAse/scripts:${PATH}"
ENV PATH="/home/apps/samtools-1.15.1:${PATH}"
ENV PATH="/home/apps/bwa-mem2-2.2.1_x64-linux:${PATH}"
ENV PATH="/home/apps/kmc/bin:${PATH}"



RUN mkdir -p /data
WORKDIR /data
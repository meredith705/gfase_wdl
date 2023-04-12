# GFAse Assembly Phasing Workflows

This repository contains WDL workflows for running [GFAse](https://github.com/rlorigro/GFAse) using proximity linkage reads or trio phasing. Several QC workflows are also provided to assess the phasing quality of the resulting haploid assemblies. 

## Inputs

Depending on which type of input files are provided to GFAse either trio-phase or proximity linkage read phasing mode is run. 
Options for input data types are: 
- HiC linked reads
- PoreC
- Parental short reads (Illumina) 
Linked read inputs are aligned to the input assembly and used to phase GFA components. Parental short reads are converted into kmers unique to each parent, matched to SNPs in the input assembly, and then used for phasing. 

## Outputs

Output files vary depending on which type of files are given as input. 
- HiC linked reads and PoreC: 
`phase0.fasta, phase1.fasta, unphased.fasta, and other files pertaining to chaining and alignment contacts.`
- Parental short reads (Illumina):
`maternal.fasta, paternal.fasta, unphased.fasta, and chaining file`

# Dependencies
- [Cromwell](https://cromwell.readthedocs.io/en/stable/)
- [jdk](https://docs.oracle.com/en/java/javase/18/install/overview-jdk-installation.html#GUID-8677A77F-231A-40F7-98B9-1FD0B48C346A)
- [docker](https://docs.docker.com/engine/install/ubuntu/)

# Installation and Usage
## Using Terra

The workflows in this repo are published to Dockstore [here](https://dockstore.org/workflows/github.com/meredith705/gfase_wdl/gfaseWorkflow) and can be easily imported to your Terra workspace. 

## Running WDL locally 

These workflows were tested locally using Cromwell which can be downloaded and installed using [these](https://cromwell.readthedocs.io/en/stable/tutorials/FiveMinuteIntro/) instructions. Here are commands to download the [current](https://github.com/broadinstitute/cromwell/releases/tag/85) release Cromwell and Womtool jar files: 
```sh
wget https://github.com/broadinstitute/cromwell/releases/download/85/cromwell-85.jar
wget https://github.com/broadinstitute/cromwell/releases/download/85/womtool-85.jar
```
A Java Development Kit ([jdk](https://docs.oracle.com/en/java/javase/18/install/overview-jdk-installation.html#GUID-8677A77F-231A-40F7-98B9-1FD0B48C346A)) is needed to run the jar files. Here is an example installation for ubuntu using apt-get:
```sh
sudo apt-get update
sudo apt install openjdk-17-jre-headless
```
Install docker ( ubuntu directions from [here](https://docs.docker.com/engine/install/ubuntu/ )) and ensure it is running:
```sh
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
  $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get update
sudo apt-get install docker-ce docker-ce-cli containerd.io
```
Clone this repo to obtain the WDL files and tiny input test files:
```sh
git clone https://github.com/meredith705/gfase_wdl.git
```
An inputs file will need to be created to run the workflow using cromwell.
For example blank input files can be created as follows:
```sh
java -jar womtool-85.jar inputs workflows/gfase.wdl > inputs.json
```
The only input documented as 'required' is the `assemblyGFA` however one of the three phasing input data types is actually required to perform any phasing. 
Input fields for each of the three input data types are as follows: 
- HiC linked reads ([json_link](https://github.com/meredith705/gfase_wdl/blob/main/inputs/input.tiny.linked_reads.json))
Paired reads are input as 2 arrays, one array for read 1 files and a second array for read 2 files. Both arrays should have the same ordering, so that read pairs are aligned together. 
```
{
  "runGFAsePhase.assemblyGFA": "tiny_test_data/test.gfa",
  "runGFAsePhase.linkedRead1Files": ["tiny_test_data/test_hic_1.fastq.gz", "tiny_test_data/test_hic_1.fastq.gz"], 
  "runGFAsePhase.linkedRead2Files": ["tiny_test_data/test_hic_2.fastq.gz", "tiny_test_data/test_hic_2.fastq.gz"],
  "runGFAsePhase.bwaAlignment.threadCount": 2
}  
```
- PoreC ([json_link](https://github.com/meredith705/gfase_wdl/blob/main/inputs/input.tiny.porec.json))
```
{
  "runGFAsePhase.assemblyGFA": "tiny_test_data/test.gfa",
  "runGFAsePhase.porecFiles": ["tiny_test_data/test_porec.fastq.gz"],
}
```
- Parental short reads ([json_link](https://github.com/meredith705/gfase_wdl/blob/main/inputs/gfase.trio.inputs.json))
```
{
  "runGFAsePhase.kmsize": 31,
  "runGFAsePhase.assemblyGFA": "tiny_test_data/test.gfa",
  "runGFAsePhase.patIlmnFiles": ["tiny_test_data/test.pat.fq"],
  "runGFAsePhase.matIlmnFiles": ["tiny_test_data/test.mat.fq"]
}
```

Cromwell is used to run the workflow using the respective inputs json file. If `CROMWELL_JAR` is the path to the cromwell jar file example commands are below.
A tiny dataset was made to help debug the WDL, it's not supposed to be test if the output is correct. 
It is intended to ensure that the WDL/commands run without an error.

```sh
## linked-read
java -jar $CROMWELL_JAR run workflows/gfase.wdl -i inputs/input.tiny.linked_reads.json

## porec
java -jar $CROMWELL_JAR run workflows/gfase.wdl -i inputs/input.tiny.porec.json

## trio
java -jar $CROMWELL_JAR run workflows/gfase.wdl -i inputs/input.tiny.trio.json

## QC - SV calling
java -jar $CROMWELL_JAR run QC/workflows/sv_evaluation.wdl -i inputs/input.tiny.qcsv.json
java -jar $CROMWELL_JAR run QC/workflows/sv_evaluation.wdl -i inputs/input.tiny.qcsv.gfa.json
java -jar $CROMWELL_JAR run QC/workflows/sv_evaluation.wdl -i inputs/input.tiny.qcsv.eval.json

## QC - QV from short reads
java -jar $CROMWELL_JAR run QC/workflows/base_qv_evaluation.wdl -i inputs/input.tiny.qcqv.json
java -jar $CROMWELL_JAR run QC/workflows/base_qv_evaluation.wdl -i inputs/input.tiny.qcqv.cram.json

## QC - QV from short reads using trio data
java -jar $CROMWELL_JAR run QC/workflows/base_qv_trio_evaluation.wdl -i inputs/input.tiny.qctrio.json

## QC - Evaluate contacts
java -jar $CROMWELL_JAR run QC/workflows/evaluate_contacts.wdl -i inputs/input.tiny.qccontacts.json
java -jar $CROMWELL_JAR run QC/workflows/evaluate_contacts.wdl -i inputs/input.tiny.qccontacts.porec.json
```

The tiny dataset was made using the python script in [tiny_test_data](tiny_test_data):

```sh
cd tiny_test_data
## simulate dataset
python3 sim_tiny_test_data.py

## gzip reads
gzip -f test_hic_1.fastq test_hic_2.fastq test_porec.fastq test.pat.fq test.mat.fq

## make a CRAM
bwa index ref.fa
bwa mem ref.fa test_hic_1.fastq.gz | samtools view -C -T ref.fa > test_reads.cram
```

# Citation
Lorig-Roach, Ryan, Melissa Meredith, Jean Monlong, Miten Jain, Hugh Olsen, Brandy McNulty, David Porubsky, et al. 2023. “Phased Nanopore Assembly with Shasta and Modular Graph Phasing with GFAse.” bioRxiv. https://doi.org/10.1101/2023.02.21.529152.



# gfase_wdl

wdl workflow to run [GFAse](https://github.com/rlorigro/GFAse) trio phasing or linked read phasing. 

Some small files to test the WDL with:

```
# gfa
wget https://rlorigro-public-files.s3.us-west-1.amazonaws.com/gfase/paolo_ul_guppy6_run14/test_subset/run14_uul_test_subset.gfa
# kmer fasta files
wget http://public.gi.ucsc.edu/~memeredith/hg002_parental_kmer_files/hg04.k31.het_hom.17-71.unique.fa
wget http://public.gi.ucsc.edu/~memeredith/hg002_parental_kmer_files/hg03.k31.het_hom_14-62.unique.fa
# linked read fastq files
wget https://rlorigro-public-files.s3.us-west-1.amazonaws.com/gfase/paolo_ul_guppy6_run14/test_subset/small_test.1.fastq
wget https://rlorigro-public-files.s3.us-west-1.amazonaws.com/gfase/paolo_ul_guppy6_run14/test_subset/small_test.2.fastq
```

Depending on input files GFAse trio-phase or linked read phasing are run. Output files also depend on whether kmer fasta files or linked read fastq files are given as input

## Testing WDL locally on tiny files

A tiny dataset was made to help debug the WDL.
It's not supposed to be test if the output is correct. 
Just that the WDL/commands run without an error.

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
```

The tiny dataset was made using the python script in [tiny_test_data](tiny_test_data):

```sh
cd tiny_test_data
## simulate dataset
python3 sim_tiny_test_data.py

## gzip reads
gzip -f test_hic_1.fastq test_hic_2.fastq test_porec.fastq

## make a CRAM
bwa index ref.fa
bwa mem ref.fa test_hic_1.fastq.gz | samtools view -C -T ref.fa > test_reads.cram
```

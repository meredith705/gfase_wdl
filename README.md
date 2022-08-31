# gfase_wdl
wdl workflow to run GFAse trio phasing or linked read phasing. 

Some small files to test the WDL with:
```
# gfa
https://rlorigro-public-files.s3.us-west-1.amazonaws.com/gfase/paolo_ul_guppy6_run14/test_subset/run14_uul_test_subset.gfa
# kmer fasta files
http://public.gi.ucsc.edu/~memeredith/hg002_parental_kmer_files/hg04.k31.het_hom.17-71.unique.fa
http://public.gi.ucsc.edu/~memeredith/hg002_parental_kmer_files/hg03.k31.het_hom_14-62.unique.fa
# linked read fastq files
https://rlorigro-public-files.s3.us-west-1.amazonaws.com/gfase/paolo_ul_guppy6_run14/test_subset/small_test.1.fastq
https://rlorigro-public-files.s3.us-west-1.amazonaws.com/gfase/paolo_ul_guppy6_run14/test_subset/small_test.2.fastq

```

Different GFAse modes are run depending on input files. Output files also depend on wether kmer fasta files or linked read fastq files are given as input

version 1.0

## KMC WDL for making kmer files from illuminia reads 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-09-20

workflow runKMC3 {
    input {
        Array[String] maternalIlmnReadFiles
        Array[String] paternalIlmnReadFiles
        Int kmerSize = 31
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 20
        Int disk_size = 4 * round(size(maternalIlmnReadFiles, 'G')) + round(size(paternalIlmnReadFiles, 'G')) + 200
        String dockerImage = "meredith705/gfase:latest"
    }

    call KMerCount { 
        input:
            maternalIlmnReadFiles = maternalIlmnReadFiles,
            paternalIlmnReadFiles = paternalIlmnReadFiles,
            memSizeGB = memSizeGB,
            threadCount = threadCount,
            disk_size = disk_size,
            dockerImage = dockerImage
    }

    output {
        File outputPaternalFa = KMerCount.paternalFasta
        File outputMaternalFa = KMerCount.maternalFasta
    }
}

task KMerCount {

    input {
        Array[String] maternalIlmnReadFiles
        Array[String] paternalIlmnReadFiles
        Int kmerSize = 31
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 20
        Int disk_size = 4 * round(size(maternalIlmnReadFiles, 'G')) + round(size(paternalIlmnReadFiles, 'G')) + 200
        String dockerImage = "meredith705/gfase:latest"
    }

    File matFiles = write_lines(maternalIlmnReadFiles)
    File patFiles = write_lines(paternalIlmnReadFiles)

    command <<<
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # make kmc tmp directory
        mkdir kmc_tmp

        # count maternal kmers
        kmc -t~{threadCount} -k~{kmerSize} @~{matFiles} maternal.kmc kmc_tmp

        # count paternal kmers
        kmc -t~{threadCount} -k~{kmerSize} @~{patFiles} paternal.kmc kmc_tmp

        # subtract kmers from each other to get unique parental kmers
        kmc_tools simple paternal.kmc maternal.kmc kmers_subtract paternal.unique.kmer
        kmc_tools simple maternal.kmc paternal.kmc kmers_subtract maternal.unique.kmer

        # paternal kmer file in to fasta files
        kmc_tools transform paternal.unique.kmer dump paternal.dump.txt
        awk 'BEGIN {OFS="_"}{msg=">pat"}{print msg,NR"\n"$1}' paternal.dump.txt > paternal.fa

        # maternal kmer file in to fasta files
        kmc_tools transform maternal.unique.kmer dump maternal.dump.txt
        awk 'BEGIN {OFS="_"}{msg=">mat"}{print msg,NR"\n"$1}' maternal.dump.txt > maternal.fa

    >>>

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File paternalFasta = "paternal.fa"
        File maternalFasta = "maternal.fa"

    }

    
}

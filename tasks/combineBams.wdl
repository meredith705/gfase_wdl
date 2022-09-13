version 1.0

## combines an array of input bams 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

task combineBams {
    input{
        Array[File] bamFiles
        # runtime configurations
        Int memSizeGB=128
        Int threadCount=64
        Int disk_size = 4 * round(size(bamFiles, 'G')) + 100
        String dockerImage="meredith705/gfase:latest"
    }
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

        samtools merge -n -@ ~{threadCount} ~{sep=" " bamFiles} -o "allBams.bam"

    >>>

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        cpu: threadCount
    }

    output {
        File outputAllbamFile = "allBams.bam"
    }
}

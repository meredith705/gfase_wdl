version 1.0

task extractAmbFromGFAse {

    input {
        File? phase0
        File? phase1
        File? unphased
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="meredith705/gfase:latest"
        Int disk_size = 5 * round(size(phase0, 'G') + size(phase1, 'G') + size(unphased, 'G')) + 20
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

        cp ~{phase0} asm0.fa
        cp ~{phase1} asm1.fa

        if [ -s ~{unphased} ]
        then
            cat ~{unphased} >> asm0.fa
            cat ~{unphased} >> asm1.fa
        fi
    >>>

    output {
        File asm0 = "asm0.fa"
        File asm1 = "asm1.fa"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

task extractAmbFromGFA {

    input {
        File? gfa
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="meredith705/gfase:latest"
        Int disk_size = 5 * round(size(gfa, 'G')) + 20
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

        python3 /home/apps/GFAse/scripts/get_haplotypes_from_shasta.py -i ~{gfa} -o haps | gzip > log.txt.gz
    >>>
    
    output {
        File asm0 = "haps/hap_0.fasta"
        File asm1 = "haps/hap_1.fasta"
        File log = "log.txt.gz"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}



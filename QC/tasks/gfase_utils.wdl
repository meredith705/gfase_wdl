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

task convertCRAMtoFASTQ {
    input{
        File cram_file
        File reference_fa
        Int memSizeGB=4
        Int threadCount=8
        Int diskSizeGB= 5 * round(size(cram_file, 'G')) + 50
        String dockerImage="quay.io/biocontainers/samtools:1.16.1--h6899075_1"
    }
    String outprefix = basename(cram_file, ".cram")
    command <<<
        set -eux -o xtrace
        
        ln -s ~{reference_fa}
        samtools fastq -@~{threadCount} --reference `basename ~{reference_fa}` -o ~{outprefix}.fq.gz ~{cram_file}
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File fastq_file = "~{outprefix}.fq.gz"
    }
}

task yakCount {
    input{
        Array[File] readFiles
        String sampleName
        Int bloomSize=37
        # runtime configurations
        Int memSizeGB=128
        Int threadCount=16
        Int diskSizeGB= 10 * round(size(readFiles, 'G')) + 50
        String dockerImage="juklucas/hpp_yak:latest"
    }
    command <<<
        set -eux -o xtrace

        # Kmer counting with https://github.com/lh3/yak.
        yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(zcat ~{sep=" " readFiles}) <(zcat ~{sep=" " readFiles})
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File yakCount = "~{sampleName}.yak"
    }
}

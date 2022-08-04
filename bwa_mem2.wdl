version 1.0

## bwa-mem2 aligner WDL for gfase proximity linked read alignment 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.ed
## 2022-08-01


task runBwaAlignment {

    input {
		File assembly_fasta
        File linked_read_fasta_1
		File linked_read_fasta_2
		# runtime configurations
		Int memSizeGB = 128
		Int threadCount = 46
		Int diskSizeGB = 64
        ## Probaly making a new docker
		String dockerImage = "meredith705/gfase:latest"
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

        # index, align, and sort reads to assembly
        bwa-mem2 index ~{assembly_fasta} \
        && \
        bwa-mem2 mem -5 -S -P ~{assembly_fasta} \
        ~{linked_read_fasta_1} \
        ~{linked_read_fasta_2} \
        | samtools sort -n -@ 24 - -o ~{assembly_fasta}_coverage.bam \

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File out_hic_to_asm_bam = ~{assembly_fasta}_coverage.bam

    }
}

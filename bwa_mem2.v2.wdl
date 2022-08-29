version 1.0

## bwa-mem2 aligner WDL for gfase proximity linked read alignment 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.ed
## 2022-08-01


workflow runBwaAlignment {
    
    input {
        File assembly_gfa_file
        Array[File?] linkedRead1Files
        Array[File?] linkedRead2Files
        File? linked_read_fasta_1_file
        File? linked_read_fasta_2_file 
        String dockerImage = "meredith705/gfase:latest" 
    }

    # run alignment
    scatter (file1 in linkedRead1Files){
        scatter (file2 in linkedRead2Files){
            call bwaAlignment {
                input:
                    assembly_gfa=assembly_gfa_file,
                    linked_read_fasta_1=file1,
                    linked_read_fasta_2=file2,
                    dockerImage=dockerImage
            }
        }
    }

    output {
        Array[File?] outputBams = flatten(bwaAlignment.outBam)

    }

}

task bwaAlignment {

    input {
		File assembly_gfa
        File? linked_read_fasta_1
		File? linked_read_fasta_2
		# runtime configurations
		Int memSizeGB = 128
		Int threadCount = 46
		Int diskSizeGB = 128
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

        if [[ ! -f "~{linked_read_fasta_1}" ]] ; then
            echo "" 

        else 
            # turn the gfa into a fasta for alignment
            python3 /home/apps/GFAse/scripts/gfa_to_fasta.py -i ~{assembly_gfa}

            # isolate the new name of the assembly fasta
            asm_fa=$(echo ~{assembly_gfa} | cut -f 1 -d ".") 
            asm_fa=$asm_fa.fasta

            # index, align, and sort reads to assembly
            bwa-mem2 index $asm_fa && \
            bwa-mem2 mem -5 -S -P $asm_fa \
            ~{linked_read_fasta_1} \
            ~{linked_read_fasta_2} \
            | samtools sort -n -@ 24 - -o "assembly.gfa.fa.coverage.bam" \

        fi

    >>>

    output {
        File? outBam = "assembly.gfa.fa.coverage.bam"

    }

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    
}

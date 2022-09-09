version 1.0

## bwa-mem2 aligner WDL for gfase proximity linked read alignment 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01


workflow runBwaAlignment {
    
    input {
        File assembly_gfa_file
        Array[File?] linkedRead1Files
        Array[File?] linkedRead2Files
        String dockerImage = "meredith705/gfase:latest" 
    }

    # run alignment
    # zip read1 and read2 into pairs to align together
    Array[Pair[File?,File?]] linked_read_pair_file_list = zip(linkedRead1Files,linkedRead2Files)
    scatter (file_pairs in linked_read_pair_file_list){
        call bwaAlignment {
            input:
                assembly_gfa=assembly_gfa_file,
                linked_read_fasta_1=file_pairs.left,
                linked_read_fasta_2=file_pairs.right,
                dockerImage=dockerImage
        }
        
    }

    output {
        Array[File?] outputBams = bwaAlignment.outBam

    }

}

task bwaAlignment {

    input {
        File assembly_gfa
        File? linked_read_fasta_1
        File? linked_read_fasta_2
        # runtime configurations
        Int memSizeGB = 185
        Int threadCount = 48
        Int disk_size = 4 * round(size(linked_read_fasta_1, 'G') + size(linked_read_fasta_2, 'G')) + 20
        String dockerImage = "meredith705/gfase:latest"
    }

    String assembly_prefix = basename(assembly_gfa, ".gfa")
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
            #python3 /home/apps/GFAse/scripts/gfa_to_fasta.py -i ~{assembly_gfa} -o assembly.fasta

            awk '/^S/{print ">"$2"\n"$3}' ~{assembly_gfa} | fold > assembly.fasta

            # store the name of the assembly fasta
            #ASM_FA=$(echo ~{assembly_gfa} | awk -F'/' '{print $(NF)}' - | cut -f 1 -d ".") 
            #ASM_FA=basename(~{assembly_gfa},'.gfa')
            #ASM_FA=$ASM_FA.fasta
            
            # store the assembly name
            ASM_FA_NAME=$(echo ~{assembly_gfa} | awk -F'/' '{print $(NF)}' | cut -f 1 -d "." )
            # store read file name
            READ1_NAME=$(echo ~{linked_read_fasta_1} | awk -F'/' '{print $(NF)}') 

            # index, align, and sort reads to assembly
            bwa-mem2 index assembly.fasta && \
            bwa-mem2 mem -t ~{threadCount} -5 -S -P assembly.fasta \
            ~{linked_read_fasta_1} \
            ~{linked_read_fasta_2} \
            | samtools sort -n -@ 24 - -o ${ASM_FA_NAME}.${READ1_NAME}.bam 

        fi

    >>>

    output {
        File? outBam = glob("*bam")[0]

    }

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    
}

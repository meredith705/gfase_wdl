version 1.0

## bwa-mem2 aligner WDL for gfase proximity linked read alignment 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01


workflow runBwaAlignment {
    
    input {
        File assembly_gfa_file
        Array[File] linkedRead1Files
        Array[File] linkedRead2Files
        String dockerImage = "meredith705/gfase:latest" 
    }

    # run alignment
    # zip read1 and read2 into pairs to align together
    Array[Pair[File,File]] linked_read_pair_file_list = zip(linkedRead1Files,linkedRead2Files)
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
        Array[File] outputBams = bwaAlignment.outBam

    }

}

task bwaAlignment {

    input {
        File assembly_gfa
        File linked_read_fasta_1
        File linked_read_fasta_2
        Int min_mapq = 1
        # runtime configurations
        Int memSizeGB = 185
        Int threadCount = 64
        Int disk_size = 15 * round(size(linked_read_fasta_1, 'G') + size(linked_read_fasta_2, 'G')) + 50
        String dockerImage = "meredith705/gfase:latest"
    }

    Int threadBwa = if threadCount < 16 then ceil(threadCount/2) else threadCount - 8
    Int threadSort = if threadCount < 16 then floor(threadCount/2) else 8
    String asm_name = sub(sub(basename(assembly_gfa), '\\.gz$', ''), '\\.gfa$', '')
    String read_name = sub(sub(sub(basename(linked_read_fasta_1), '\\.gz$', ''), '\\.fq$', ''), '\\.fastq$', '')
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

        # turn the gfa into a fasta for alignment
        python3 /home/apps/GFAse/scripts/gfa_to_fasta.py -i ~{assembly_gfa} -o ./assembly.fasta

        # index, align, and sort reads to assembly
        bwa-mem2 index assembly.fasta
        bwa-mem2 mem -t ~{threadBwa} -5 -S -P assembly.fasta \
                 ~{linked_read_fasta_1} ~{linked_read_fasta_2} | \
            samtools view -h -q ~{min_mapq} - | \
            samtools sort -n -@ ~{threadSort} - -O BAM -o ~{asm_name}.~{read_name}.bam 
    >>>

    output {
        File outBam = "${asm_name}.${read_name}.bam "
    }

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    
}

task diploidBwaAlignment {

    input {
        File haps_fasta
        File linked_read_fasta_1
        File linked_read_fasta_2
        Int min_mapq = 1
        # runtime configurations
        Int memSizeGB = 185
        Int threadCount = 64
        Int disk_size = 15 * round(size(linked_read_fasta_1, 'G') + size(linked_read_fasta_2, 'G')) + 50
        String dockerImage = "meredith705/gfase:latest"
    }

    Int threadBwa = if threadCount < 16 then ceil(threadCount/2) else threadCount - 8
    Int threadSort = if threadCount < 16 then floor(threadCount/2) else 8
    String ref_name = sub(sub(basename(haps_fasta), '\\.fasta$', ''), '\\.fa$', '')
    String read_name = sub(sub(sub(basename(linked_read_fasta_1), '\\.gz$', ''), '\\.fq$', ''), '\\.fastq$', '')
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

        ln -s ~{haps_fasta} ref.fasta
        
        # index, align, and sort reads to assembly
        bwa-mem2 index ref.fasta
        bwa-mem2 mem -t ~{threadBwa} -5 -S -P ref.fasta \
                 ~{linked_read_fasta_1} ~{linked_read_fasta_2} | \
            samtools view -h -q ~{min_mapq} - | \
            samtools sort -n -@ ~{threadSort} - -O BAM -o ~{ref_name}.~{read_name}.bam 
    >>>

    output {
        File outBam = "${ref_name}.${read_name}.bam "
    }

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }    
}

task minimap2Alignment {

    input {
        File assembly_gfa
        File porec_file
        Int min_mapq = 1
        # runtime configurations
        Int memSizeGB = 185
        Int threadCount = 64
        Int disk_size = 5 * round(size(porec_file, 'G')) + 50
        String dockerImage = "meredith705/gfase:latest"
    }

    Int threadAlign = if threadCount < 16 then ceil(threadCount/2) else threadCount - 8
    Int threadView = if threadCount < 16 then floor(threadCount/2) else 8
    String asm_name = sub(sub(basename(assembly_gfa), '\\.gz$', ''), '\\.gfa$', '')
    String read_name = sub(sub(sub(basename(porec_file), '\\.gz$', ''), '\\.fq$', ''), '\\.fastq$', '')
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

        # turn the gfa into a fasta for alignment
        python3 /home/apps/GFAse/scripts/gfa_to_fasta.py -i ~{assembly_gfa} -o ./assembly.fasta

        # align with minimap2
        minimap2 -a -x map-ont -k 17 -t ~{threadAlign} \
                 -K 10g -I 8g \
                 assembly.fasta \
                 ~{porec_file} | samtools view -bh -@ ~{threadView} -q ~{min_mapq} -o ~{asm_name}.~{read_name}.bam -O BAM -
    >>>

    output {
        File outBam = "${asm_name}.${read_name}.bam "
    }

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    
}

task diploidMinimap2Alignment {

    input {
        File haps_fasta
        File porec_file
        Int min_mapq = 1
        # runtime configurations
        Int memSizeGB = 185
        Int threadCount = 64
        Int disk_size = 5 * round(size(porec_file, 'G')) + 50
        String dockerImage = "meredith705/gfase:latest"
    }

    Int threadAlign = if threadCount < 16 then ceil(threadCount/2) else threadCount - 8
    Int threadView = if threadCount < 16 then floor(threadCount/2) else 8
    String ref_name = sub(sub(basename(haps_fasta), '\\.fasta$', ''), '\\.fa$', '')
    String read_name = sub(sub(sub(basename(porec_file), '\\.gz$', ''), '\\.fq$', ''), '\\.fastq$', '')
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
      
        # align with minimap2
        minimap2 -a -x map-ont -k 17 -t ~{threadAlign} \
                 -K 10g -I 8g \
                 ~{haps_fasta} \
                 ~{porec_file} | samtools view -bh -@ ~{threadView} -q ~{min_mapq} -o ~{ref_name}.~{read_name}.bam -O BAM -
    >>>

    output {
        File outBam = "${ref_name}.${read_name}.bam "
    }

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    
}


version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

import "../tasks/bwa_mem2.wdl" as bwamem_t

workflow runGFAsePhase {
	
    input {
        File assemblyGFA                            # Input assembly GFA file from assembler
        #### trio data input ########
        File? patKmerFile                           # Input Fasta file of paternal kmers ( eg. made with KMC3 )
        File? matKmerFile                           # Input Fasta file of maternal kmers ( eg. made with KMC3 )
        Int kmsize = 31                            # Size of kmers in the input files
        #### linked read input ########
        Array[File] linkedRead1Files = []              # Input array of read1 fastq files from a set of paired linked reads
        Array[File] linkedRead2Files = []              # Input array of read2 fastq files from a set of paired linked reads 
        String dockerImage = "meredith705/gfase:latest" 
    }

    # if linked reads are supplied as input align linked reads and phase the gfa
    if (length(linkedRead1Files) > 0 && length(linkedRead2Files) > 0){

        # zip read1 and read2 into pairs to align together
        Array[Pair[File,File]] linked_read_pair_file_list = zip(select_all(linkedRead1Files),select_all(linkedRead2Files))
        scatter (file_pairs in linked_read_pair_file_list){
            call bwamem_t.bwaAlignment {
                input:
                    assembly_gfa        = assemblyGFA,
                    linked_read_fasta_1 = file_pairs.left,
                    linked_read_fasta_2 = file_pairs.right,
                    dockerImage         = dockerImage
            }
        }

        # phase the gfa using the linked read alignment
        call gfaseLinkedRead {
            input:
                assemblyGfa         = assemblyGFA,
                bamFiles            = bwaAlignment.outBam,
                dockerImage         = dockerImage
        }

            
    }

    if (defined(patKmerFile) && defined(matKmerFile)){
        # trio phase gfa
        call gfaseTrioPhase {
            input:
                assemblyGfa         = assemblyGFA,
                patKmerFa           = patKmerFile,
                matKmerFa           = matKmerFile,
                kSize               = kmsize,
                dockerImage         = dockerImage
        }

    }

    output {
        File? outputMatAssembly             = gfaseTrioPhase.outMatAssembly
        File? outputPatAssembly             = gfaseTrioPhase.outPatAssembly
        File? outputUnphasedAssembly        = gfaseTrioPhase.outUnphasedAssembly
        File? outputPhaseChains             = gfaseTrioPhase.outPhaseChains
        File? outputConfig                  = gfaseLinkedRead.outconfig
        File? outputContacts                = gfaseLinkedRead.outcontacts 
        File? outputFaPhase0                = gfaseLinkedRead.outFastaP0 
        File? outputFaPhase1                = gfaseLinkedRead.outFastaP1 
        File? outputPhases                  = gfaseLinkedRead.outPhases 
        File? outputFaUnphasedInnitial      = gfaseLinkedRead.outFastaUnP 
        File? outputFaUnPin                 = gfaseLinkedRead.outFastaUnPinn
    }

}

task gfaseTrioPhase {
    input {
        File assemblyGfa
        # trio parameters
        File? patKmerFa 
        File? matKmerFa
        Int kSize = 0
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 1
        Int disk_size = 4 * round(size(assemblyGfa, 'G')) + round(size(patKmerFa, 'G')) + round(size(matKmerFa, 'G')) + 100
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
        set -o xtrace

        echo trio phase >> log.txt
        # Phase the gfa
        phase_haplotype_paths \
        -k ~{kSize} \
        -i ~{assemblyGfa} \
        -p ~{patKmerFa} \
        -m ~{matKmerFa} \
        -o gfase
        
        echo done >> log.txt
    >>>

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outMatAssembly  = "gfase/maternal.fasta"
        File outPatAssembly  = "gfase/paternal.fasta"
        File outUnphasedAssembly = "gfase/unphased.fasta"
        File outPhaseChains  = "gfase/phase_chains.csv"
        File outWDLlog       = "log.txt"

    }
}

task gfaseLinkedRead {
    input {
        File assemblyGfa
        # hi-c phasing bam
        Array[File] bamFiles
        Int min_mapq = 2
        String target_contig_prefix = "PR"
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 46
        Int disk_size = 10 * round(size(assemblyGfa, 'G')) + round(size(bamFiles, 'G')) + 100
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

        samtools merge -n -@ ~{threadCount} -o allBams.bam ~{sep=" " bamFiles}

        echo hi-c phase >> log.txt
        # hi-c 
        phase_contacts \
        -i allBams.bam \
        -g ~{assemblyGfa} \
        -o gfase \
        -m ~{min_mapq} \
        -t ~{threadCount} \
        -p ~{target_contig_prefix}
                
        echo done >> log.txt
 
    >>>

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outconfig       = "gfase/config.csv"
        File outcontacts     = "gfase/contacts.csv"
        File outFastaP0      = "gfase/phase_0.fasta"
        File outFastaP1      = "gfase/phase_1.fasta"
        File outPhases       = "gfase/phases.csv"
        File outFastaUnP     = "gfase/unphased.fasta"
        File outFastaUnPinn  = "gfase/unphased_initial.fasta"
        File outWDLlog       = "log.txt"

    }
}




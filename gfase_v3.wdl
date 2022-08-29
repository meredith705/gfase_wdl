version 1.0

## GFAse WDL for gfa phasing, using trio or linked reads 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.ed
## 2022-08-01

import "bwa_mem2.wdl" as bwamem_t

workflow runGFAsePhase {
	
    input {
        File assemblyGFA
        File? patKmerFile
        File? matKmerFile
        #### work on this ########
        Boolean linked_reads = "True"
        Array[File?] linkedRead1Files
        Array[File?] linkedRead2Files
        File? linked_read_fasta_1_File
        File? linked_read_fasta_2_File
        Int? kmsize 
        String dockerImage = "meredith705/gfase:latest" 
    }

    # align linked reads - scattered Pairs of reads
    if (linked_reads){
        scatter (file1 in linkedRead1Files){
            scatter (file2 in linkedRead2Files){
                call bwamem_t.bwaAlignment as bwamem{
                    input:
                        assembly_gfa        = assemblyGFA,
                        linked_read_fasta_1 = file1,
                        linked_read_fasta_2 = file2,
                        dockerImage         = dockerImage
                }

            }
         
        }

        # combine the bam's 
        call combineBams as cb {
            input:
                bamFiles =  flatten(bwamem.outBam)
        }
    }

    
    


    # phase gfa
    call gfasePhase {
        input:
            assemblyGfa         = assemblyGFA,
            patKmerFa           = patKmerFile,
            matKmerFa           = matKmerFile,
            kSize               = kmsize,
            linked_reads_bool   = linked_reads,
            alignmentBam        = cb.outputAllbamFile,
            dockerImage         = dockerImage
    }

    output {
        File? outputMatAssembly = gfasePhase.outMatAssembly
        File? outputPatAssembly = gfasePhase.outPatAssembly
        File? outputUnPAssembly = gfasePhase.outUnphasedAssembly
        File? outputPhaseChains = gfasePhase.outPhaseChains
        File? outputConfig    = gfasePhase.outconfig
        File? outputContacts  = gfasePhase.outcontacts 
        File? outputFaP0      = gfasePhase.outFastaP0 
        File? outputFaP1      = gfasePhase.outFastaP1 
        File? outputPhases    = gfasePhase.outPhases 
        File? outputFaUnP     = gfasePhase.outFastaUnP 
        File? outputFaUnPin   = gfasePhase.outFastaUnPinn
    }

}

task combineBams {
    input{
        Array[File?] bamFiles
        # runtime configurations
        Int memSizeGB=128
        Int threadCount=16
        Int diskSizeGB=256
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

        samtools merge -n -@ 24 <(cat ~{sep=" " bamFiles}) -o "allBams.bam"

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File outputAllbamFile = "allBams.bam"
    }
}

task gfasePhase {
    input {
        File assemblyGfa
        # trio parameters
        File? patKmerFa 
        File? matKmerFa
        Int? kSize = 0
        # hi-c phasing bam
        Boolean linked_reads_bool 
        File? alignmentBam
        Int m = 2
        String p = "PR"
        Int t = 46
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 1
        Int diskSizeGB = 64
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

        echo checking for alignmentBam. >> log.txt
        if [[ ! "~{linked_reads_bool}" ]] ; then
            echo no alignmentBam checking for parental kmer files. >> log.txt
            if [[ -f "~{patKmerFa}" ]] ; then
                echo trio phase >> log.txt
                # Phase the gfa
                phase_haplotype_paths \
                -k ~{kSize} \
                -i ~{assemblyGfa} \
                -p ~{patKmerFa} \
                -m ~{matKmerFa}
            else
                echo no trio files. >> log.txt
                echo no phasing can be done without either trio kmer fasta files or proximity liagation fasta files. >> log.txt
            fi
            
        else
            echo hi-c phase >> log.txt
            # hi-c 
            phase_contacts \
            -i ~{alignmentBam} \
            -g ~{assemblyGfa} \
            -o gfase \
            -m ~{m} \
            -p ~{p} 
            #output:
            #config.csv  contacts.csv  phase_0.fasta  phase_1.fasta  phase_chains.csv  phases.csv  unphased.fasta  unphased_initial.fasta

            tar czvf gfase.tar.gz gfase

            
        fi

        
        echo done >> log.txt
 

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
    }

    output {
        File? outMatAssembly  = "maternal.fasta"
        File? outPatAssembly  = "paternal.fasta"
        File? outUnphasedAssembly = "unphased.fasta"
        File? outPhaseChains  = "phase_chains.csv"
        File? outconfig       = "gfase/config.csv"
        File? outcontacts     = "gfase/contacts.csv"
        File? outFastaP0      = "gfase/phase_0.fasta"
        File? outFastaP1      = "gfase/phase_1.fasta"
        File? outPhases       = "gfase/phases.csv"
        File? outFastaUnP     = "gfase/unphased.fasta"
        File? outFastaUnPinn  = "gfase/unphased_initial.fasta"
        File? outWDLlog       = "log.txt"

    }
}




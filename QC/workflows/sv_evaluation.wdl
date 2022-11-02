version 1.0

import "https://raw.githubusercontent.com/nanoporegenomics/variant_calling_wf/main/wdl/tasks/dipdiff.wdl" as dipdiff_t

workflow gfaseSVevaluation {

    input {
        File? assemblyGfa
        File? assemblyPhase0
        File? assemblyPhase1
        File? assemblyUnphased
        File? truthVcf
        File reference
        File vntrAnnotations
        String sampleId="sample"
    }

    if(defined(assemblyGfa)){
        call extractAmbFromGFA {
            input:
            gfa=assemblyGfa
        }
    }
    if(!defined(assemblyGfa) && defined(assemblyPhase1) && defined(assemblyPhase0)){
        call extractAmbFromGFAse {
            input:
            phase0=assemblyPhase0,
            phase1=assemblyPhase1,
            unphased=assemblyUnphased
        }
    }
    
    File asm0 = select_first([extractAmbFromGFA.asm0, extractAmbFromGFAse.asm0])
    File asm1 = select_first([extractAmbFromGFA.asm1, extractAmbFromGFAse.asm1])
    
    call dipdiff_t.dipdiff_t {
        input:
        ctgsPat=asm0,
        ctgsMat=asm1,
        reference=reference,
        vntrAnnotations=vntrAnnotations
    }

    # if(defined(truthVcf)){
    #     # evaluation with truvari?
    # }
    
    output {
		File vcfUnphased = dipdiff_t.dipdiffUnphasedVcf
		File vcfPhased = dipdiff_t.dipdiffPhasedVcf
    }
}


task extractAmbFromGFAse {

    input {
        File? phase0
        File? phase1
        File? unphased
        Int threadCount = 1
        Int memoryGB = 4
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

        cp ~{phase0} asm0.fa
        cp ~{phase1} asm1.fa

        if [ -s ~{unphased} ]
        then
            cat ~{unphased} >> asm0.fa
        fi
    >>>

    output {
        File asm0 = "asm0.fa"
        File asm1 = "asm1.fa"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: memoryGB + "GB"
    }
}

task extractAmbFromGFA {

    input {
        File? gfa
        Int threadCount = 1
        Int memoryGB = 4
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

        python3 /home/apps/GFAse/scripts/get_haplotypes_from_shasta.py -i ~{gfa} -o haps
    >>>
    
    output {
        File asm0 = "haps/hap_0.fasta"
        File asm1 = "haps/hap_1.fasta"
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: memoryGB + "GB"
    }
}



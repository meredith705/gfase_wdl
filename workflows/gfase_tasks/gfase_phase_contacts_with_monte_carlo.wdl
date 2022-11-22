version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

workflow RunGfaseLinkedReadmc {
    input {
        File assemblyGFA                        # gfa from Shasta/verrko
        Array[File] bamFiles                    # hi-c phasing bam
        Int min_mapq = 1                        # default min mapq is 1
        Int threadCount = 46                    # Minimum required mapq value for mapping to be counted.
        String otherPhaseContactsArugments = "" # for verrko graphs use --use_homology --skip_unzip

        # runtime configurations
        Int memSizeGB = 128
        Int disk_size = 10 * round(size(assemblyGFA, 'G') + size(bamFiles, 'G')) + 100
        String dockerImage = "meredith705/gfase:latest"
    }

    # phase the gfa using the linked read alignment
    call gfase_phase_contacts_with_monte_carlo as gfasePhaseContactsMC {
        input:
            assemblyGfa         = assemblyGFA,
            bamFiles            = bamFiles,
            dockerImage         = dockerImage,
            otherPhaseContactsArugments = otherPhaseContactsArugments
        }

    output {
        File outputConfig                  = gfasePhaseContactsMC.outconfig
        File outputContacts                = gfasePhaseContactsMC.outcontacts 
        File outputFaPhase0                = gfasePhaseContactsMC.outFastaP0 
        File outputFaPhase1                = gfasePhaseContactsMC.outFastaP1 
        File outputPhases                  = gfasePhaseContactsMC.outPhases 
        File? outputUnzipGFA                = gfasePhaseContactsMC.outUnzipGFA
        File outputChainedGFA              = gfasePhaseContactsMC.outChainedGFA
        File outputFaUnphased              = gfasePhaseContactsMC.outFastaUnP 
    }

}

task gfase_phase_contacts_with_monte_carlo {
    input {
        File assemblyGfa                     # gfa from Shasta/verrko
        Array[File] bamFiles                 # hi-c phasing bam
        Int min_mapq = 1                     # default min mapq is 1
        Int threadCount = 46                 # Minimum required mapq value for mapping to be counted.
        String otherPhaseContactsArugments = ""

        # other possible arugments
        Int? core_iterations                 # Iterations for each shallow convergence in the sampling process (uses 3*core_iterations)
        Int? sample_size                     # Num shallowly converged phase states to sample from
        Int? n_rounds                        # Rounds to sample and merge

        # runtime configurations
        Int memSizeGB = 128
        Int disk_size = 10 * round(size(assemblyGfa, 'G') + size(bamFiles, 'G')) + 100
        String dockerImage = "meredith705/gfase:latest"
    }

    command <<<
        set -o pipefail
        set -e
        set -u
        set -o xtrace

        samtools merge -n -@ ~{threadCount} -o allBams.bam ~{sep=" " bamFiles}

        echo hi-c phase >> log.txt
        # hi-c 
        phase_contacts_with_monte_carlo \
        -i allBams.bam \
        -g ~{assemblyGfa} \
        -o gfase \
        -m ~{min_mapq} \
        -t ~{threadCount} \
        ~{otherPhaseContactsArugments}
                
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
        # unzipped gfa and chained.gfa
        File? outUnzipGFA     = "gfase/unzipped.gfa"
        File outChainedGFA   = "gfase/chained.gfa"
        File outFastaUnP     = "gfase/unphased.fasta"
        File outWDLlog       = "log.txt"
        # removed
        #File outFastaUnPinn  = "gfase/unphased_initial.fasta"
        

    }
}




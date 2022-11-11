version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

workflow RunGfaseLinkedReadmc {
    input {
        File assemblyGFA_w                    # gfa from Shasta/verrko
        Array[File] bamFiles_w                 # hi-c phasing bam
        Int min_mapq_w = 1                     # default min mapq is 1
        Int threadCount_w = 46                 # Minimum required mapq value for mapping to be counted.
        String otherPhaseContactsArugments_w = "--use_homology --skip_unzip"

        # other possible arugments
        Int? core_iterations_w                 # Iterations for each shallow convergence in the sampling process (uses 3*core_iterations)
        Int? sample_size_w                     # Num shallowly converged phase states to sample from
        Int? n_rounds_w                        # Rounds to sample and merge

        # runtime configurations
        Int memSizeGB_w = 128
        Int disk_size_w = 10 * round(size(assemblyGfa, 'G')) + round(size(bamFiles, 'G')) + 100
        String dockerImage_w = "meredith705/gfase:latest"
    }

    # phase the gfa using the linked read alignment
    call gfase_phase_contacts_with_monte_carlo as gfaseLinkedRead {
        input:
            assemblyGfa         = assemblyGFA_w,
            bamFiles            = bamFiles_w,
            dockerImage         = dockerImage_w
        }

            
    }


    output {
        File outputConfig                  = gfaseLinkedRead.outconfig
        File outputContacts                = gfaseLinkedRead.outcontacts 
        File outputFaPhase0                = gfaseLinkedRead.outFastaP0 
        File outputFaPhase1                = gfaseLinkedRead.outFastaP1 
        File outputPhases                  = gfaseLinkedRead.outPhases 
        File outputUnzipGFA                = gfaseLinkedRead.outUnzipGFA
        File outputChainedGFA              = gfaseLinkedRead.outChainedGFA
        File outputFaUnphased              = gfaseLinkedRead.outFastaUnP 
    }

}

task gfase_phase_contacts_with_monte_carlo {
    input {
        File assemblyGfa                    # gfa from Shasta/verrko
        Array[File] bamFiles                 # hi-c phasing bam
        Int min_mapq = 1                     # default min mapq is 1
        Int threadCount = 46                 # Minimum required mapq value for mapping to be counted.
        String otherPhaseContactsArugments = "--use_homology --skip_unzip"

        # other possible arugments
        Int? core_iterations                 # Iterations for each shallow convergence in the sampling process (uses 3*core_iterations)
        Int? sample_size                     # Num shallowly converged phase states to sample from
        Int? n_rounds                        # Rounds to sample and merge

        # runtime configurations
        Int memSizeGB = 128
        Int disk_size = 10 * round(size(assemblyGfa, 'G')) + round(size(bamFiles, 'G')) + 100
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
        -t ~{threadCount} 
                
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
        File outUnzipGFA     = "gfase/unzipped.gfa"
        File outChainedGFA   = "gfase/chained.gfa"
        File outFastaUnP     = "gfase/unphased.fasta"
        File outWDLlog       = "log.txt"
        # removed
        #File outFastaUnPinn  = "gfase/unphased_initial.fasta"
        

    }
}




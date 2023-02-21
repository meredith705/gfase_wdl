version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

import "../../tasks/kmc.wdl" as kmc_t

workflow RunGFAseTrioPhase {
    input {
        File assemblyGfa                         # gfa from Shasta/verrko
        Array[File] patKmerFa                    # array of paternal files
        Array[File] matKmerFa                    # array of maternal files
        Int kSize = 31                           # kmer size for parents and child

        # runtime configurations
        Int memSizeGB = 128
        Int disk_size = 10 * round(size(assemblyGfa, 'G')) + 2 * round(size(patKmerFa, 'G')) + 2 * round(size(matKmerFa, 'G')) + 100
        String dockerImage = "meredith705/gfase:latest"
    }


    # turn the fastas into kmers for phasing 
    call kmc_t.KMerCount {
        input:
            maternalIlmnReadFiles = matKmerFa,
            paternalIlmnReadFiles = patKmerFa,
            kmerSize              = kSize
    }

    # phase the gfa using the parental kmers
    call gfaseTrioPhase {
        input:
            assemblyGfa         = assemblyGfa,
            patKmers           = KMerCount.paternalFasta,
            matKmers           = KMerCount.maternalFasta,
            kSize               = kSize,
            dockerImage         = dockerImage
        }

    output {
        File? outputMatAssembly             = gfaseTrioPhase.outMatAssembly
        File? outputPatAssembly             = gfaseTrioPhase.outPatAssembly
        File? outputTrioUnphasedAssembly    = gfaseTrioPhase.outUnphasedAssembly
        File? outputPhaseChains             = gfaseTrioPhase.outPhaseChains
    }

}


task gfaseTrioPhase {
    input {
        File assemblyGfa
        # trio parameters
        File patKmers 
        File matKmers
        Int kSize = 0
        # runtime configurations
        Int memSizeGB = 128
        Int threadCount = 1
        Int disk_size = 4 * round(size(assemblyGfa, 'G')) + round(size(patKmers, 'G')) + round(size(matKmers, 'G')) + 100
        String dockerImage = "meredith705/gfase:latest"
        Int preemptible = 2
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
        -p ~{patKmers} \
        -m ~{matKmers} \
        -o gfase
        
        echo done >> log.txt
    >>>

    runtime {
        docker: dockerImage
        disks: "local-disk " + disk_size + " SSD"
        memory: memSizeGB + " GB"
        cpu: threadCount
        preemptible: preemptible
    }

    output {
        File outMatAssembly      = "gfase/maternal.fasta"
        File outPatAssembly      = "gfase/paternal.fasta"
        File outUnphasedAssembly = "gfase/unphased.fasta"
        File outPhaseChains      = "gfase/phase_chains.csv"
        File outWDLlog           = "log.txt"

    }
}
version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

workflow RunGFAseTrioPhase {
    input {
        File assemblyGfa                         # gfa from Shasta/verrko
        File? patKmerFa                          # paternal file
        File? matKmerFa                          # maternal file
        Int kSize = 31                           # kmer size for parents and child

        # runtime configurations
        Int memSizeGB = 128
        Int disk_size = 10 * round(size(assemblyGfa, 'G')) + round(size(patKmerFa, 'G')) + + round(size(matKmerFa, 'G')) + 100
        String dockerImage = "meredith705/gfase:latest"
    }

    # phase the gfa using the linked read alignment
    call gfaseTrioPhase {
        input:
            assemblyGfa         = assemblyGfa,
            patKmerFa           = patKmerFa,
            matKmerFa           = matKmerFa,
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
        File outMatAssembly      = "gfase/maternal.fasta"
        File outPatAssembly      = "gfase/paternal.fasta"
        File outUnphasedAssembly = "gfase/unphased.fasta"
        File outPhaseChains      = "gfase/phase_chains.csv"
        File outWDLlog           = "log.txt"

    }
}
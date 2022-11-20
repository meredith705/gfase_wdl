version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

import "../tasks/proximityAlign.wdl" as align_t
import "gfase_tasks/gfase_phase_contacts_with_monte_carlo.wdl" as gfaseLR

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
        #### porec input
        Array[File] porecFiles = []                           # Input array of fastq files from poreC experiments
        String dockerImage = "meredith705/gfase:latest" 
    }

    # if proximity data is supplied as input align reads and phase the gfa
    if (length(porecFiles) > 0 || (length(linkedRead1Files) > 0 && length(linkedRead2Files) > 0)){

        # if proximity data is linked reads from HiC
        if (length(linkedRead1Files) > 0 && length(linkedRead2Files) > 0){
            # zip read1 and read2 into pairs to align together
            Array[Pair[File,File]] linked_read_pair_file_list = zip(select_all(linkedRead1Files),select_all(linkedRead2Files))
            scatter (file_pairs in linked_read_pair_file_list){
                call align_t.bwaAlignment {
                    input:
                    assembly_gfa        = assemblyGFA,
                    linked_read_fasta_1 = file_pairs.left,
                    linked_read_fasta_2 = file_pairs.right,
                    dockerImage         = dockerImage
                }
            }
        }

        # if proximity data is porec reads
        if (length(porecFiles) > 0){
            # zip read1 and read2 into pairs to align together
            scatter (porecF in porecFiles){
                call align_t.minimap2Alignment {
                    input:
                    assembly_gfa = assemblyGFA,
                    porec_file = porecF
                }
            }
        }
        
        Array[File] proxBam = select_first([bwaAlignment.outBam, minimap2Alignment.outBam])

        # phase the gfa using the linked read alignment
        call gfaseLR.RunGfaseLinkedReadmc as gfaseLinkedRead{
            input:
            assemblyGFA         = assemblyGFA,
            bamFiles            = proxBam,
            dockerImage         = dockerImage
        }        
    }
    
    # otherwise trio-phasing based on kmers
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
        File? outputTrioUnphasedAssembly    = gfaseTrioPhase.outUnphasedAssembly
        File? outputPhaseChains             = gfaseTrioPhase.outPhaseChains
        File? outputConfig                  = gfaseLinkedRead.outputConfig
        File? outputContacts                = gfaseLinkedRead.outputContacts
        File? outputFaPhase0                = gfaseLinkedRead.outputFaPhase0 
        File? outputFaPhase1                = gfaseLinkedRead.outputFaPhase1 
        File? outputPhases                  = gfaseLinkedRead.outputPhases 
        File? outputUnzipGFA                = gfaseLinkedRead.outputUnzipGFA
        File? outputChainedGFA              = gfaseLinkedRead.outputChainedGFA
        File? outputLRFaUnphased            = gfaseLinkedRead.outputFaUnphased
        
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


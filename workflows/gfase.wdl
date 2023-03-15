version 1.0

## gfase wdl 
## Maintainer: Melissa Meredith 
## mmmeredi@ucsc.edu
## 2022-08-01

import "../tasks/proximityAlign.wdl" as align_t
import "gfase_tasks/gfase_phase_contacts_with_monte_carlo.wdl" as gfaseLR
import "gfase_tasks/gfase_trio_phase.wdl" as gfaseTrio

workflow runGFAsePhase {
	
    input {
        File assemblyGFA                            # Input assembly GFA file from assembler
        #### trio data input ########
        Array[File] patIlmnFiles = []               # Input array of parental fasta files 
        Array[File] matIlmnFiles = []               # Input array of marental fasta files 
        Int kmsize = 31                             # Size of kmers in the input files
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
    if (length(patIlmnFiles) > 0 && length(matIlmnFiles) > 0){
        # trio phase gfa
        call gfaseTrio.RunGFAseTrioPhase as gfaseTrioPhase {
            input:
                assemblyGfa         = assemblyGFA,
                patKmerFa           = patIlmnFiles,
                matKmerFa           = matIlmnFiles,
                kSize               = kmsize,
                dockerImage         = dockerImage
        }

    }

    output {
        File? outputMatAssembly             = gfaseTrioPhase.outputMatAssembly
        File? outputPatAssembly             = gfaseTrioPhase.outputPatAssembly
        File? outputTrioUnphasedAssembly    = gfaseTrioPhase.outputTrioUnphasedAssembly
        File? outputPhaseChains             = gfaseTrioPhase.outputPhaseChains
        File? outputConfig                  = gfaseLinkedRead.outputConfig
        File? outputContacts                = gfaseLinkedRead.outputContacts
        File? outputFaPhase0                = gfaseLinkedRead.outputFaPhase0 
        File? outputFaPhase1                = gfaseLinkedRead.outputFaPhase1 
        File? outputPhases                  = gfaseLinkedRead.outputPhases 
        File? outputUnzipGFA                = gfaseLinkedRead.outputUnzipGFA
        File? outputChainedGFA              = gfaseLinkedRead.outputChainedGFA
        File? outputLRFaUnphased            = gfaseLinkedRead.outputFaUnphased
        Array[File]? outputLinkedReadBam    = proxBam
    }

}




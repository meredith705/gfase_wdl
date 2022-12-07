version 1.0

import "../../tasks/proximityAlign.wdl" as align_t


workflow evaluateContacts {

    meta {
	    author: "Jean Monlong"
        email: "jmonlong@ucsc.edu"
        description: "Compute some stats about contacts in a phased assembly from proximity sequencing (HiC or PoreC)"
    }

    parameter_meta {
        linkedReads1: "First pairs of HiC reads (gzipped FASTQ). Optional. Use with refHap1/2. Use one of linkedReads1/2 or porecReads."
        linkedReads2: "Second pairs of HiC reads (gzipped FASTQ). Optional. Use with refHap1/2. Use one of linkedReads1/2 or porecReads."
        porecReads: "Long reads, e.g. PoreC (gzipped FASTQ). Optional. Use with refHap1/2. Use one of linkedReads1/2 or porecReads."
        refHap1: "Haplotype 1 FASTA. Optional. Use with linkedReads1/2 or porecReads."
        refHap2: "Haplotype 2 FASTA. Optional. Use with linkedReads1/2 or porecReads."
        bamFile: "Reads aligned to the phased assembly (BAM file). Optional. Use with phaseContactsCsv, and instead of linkedReads1/2/porecReads and refHap1/2"
        phaseContactsCsv: "CSV file defining the assigned phases of contigs in the bamFile. Optional. Use with bamFile, and instead of linkedReads1/2/porecReads and refHap1/2"
    }

    input {
        File? linkedReads1
        File? linkedReads2
        File? porecReads
        File? refHap1
        File? refHap2
        File? bamFile
        File? phaseContactsCsv
    }

    if (!defined(phaseContactsCsv) && !defined(bamFile)){
        call mergePhasedFasta {
            input:
            hap1_fasta=refHap1,
            hap2_fasta=refHap2
        }
        
        # if proximity data is linked reads from HiC
        if (defined(linkedReads1) && defined(linkedReads2)){
            call align_t.diploidBwaAlignment {
                input:
                linked_read_fasta_1 = select_first([linkedReads1]),
                linked_read_fasta_2 = select_first([linkedReads2]),
                haps_fasta=mergePhasedFasta.outMergedHaps
            }
        }
    
        # if proximity data is porec reads
        if (defined(porecReads)){
            call align_t.diploidMinimap2Alignment {
                input:
                porec_file = select_first([porecReads]),
                haps_fasta=mergePhasedFasta.outMergedHaps
            }
        }
    }

    File contacts = select_first([phaseContactsCsv, mergePhasedFasta.outContacts])
    File proxBam = select_first([bamFile, diploidBwaAlignment.outBam, diploidMinimap2Alignment.outBam])
    
    # phase the gfa using the linked read alignment
    call evaluateContactsTask {
        input:
        bam=proxBam,
        contacts=contacts
    }
       
    output {
		File outputTar = evaluateContactsTask.outputTar
		File contactsSummary = evaluateContactsTask.contactsSummary
		File phasingSummary = evaluateContactsTask.phasingSummary
    }
}

##
#### TASKS
##

task mergePhasedFasta {
    input {
        File? hap1_fasta
        File? hap2_fasta
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="meredith705/gfase:latest"
        Int disk_size = 3 * round(size(hap1_fasta, 'G') + size(hap2_fasta, 'G')) + 20
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

        python3 /home/apps/GFAse/scripts/merge_phased_fasta.py -a ~{hap1_fasta} -b  ~{hap2_fasta} -o merged_haps/
        
    >>>
    output {
        File outContacts = "merged_haps/phases.csv"
        File outMergedHaps = "merged_haps/merged.fasta"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}


task evaluateContactsTask {
    input {
        File bam
        File contacts
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="meredith705/gfase:latest"
        Int disk_size = 3 * round(size(bam, 'G')) + 20
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

        evaluate_contacts -i ~{bam} -p ~{contacts} -o contact_eval_results

        tar -czvf contact_eval_results.tar.gz contact_eval_results/*
    >>>
    output {
        File outputTar = "contact_eval_results.tar.gz"
        File contactsSummary = "contact_eval_results/contacts_summary.csv"
        File phasingSummary = "contact_eval_results/phasing_summary.csv"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

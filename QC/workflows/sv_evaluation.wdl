version 1.0

import "../tasks/gfase_utils.wdl" as utils
import "https://raw.githubusercontent.com/nanoporegenomics/variant_calling_wf/main/wdl/tasks/dipdiff.wdl" as dipdiff_t

workflow gfaseSVevaluation {

    input {
        File? assemblyGfa
        File? assemblyPhase0
        File? assemblyPhase1
        File? assemblyUnphased
        File? callVcf
        File? truthVcf
        File? confidentRegions
        File reference
        File? vntrAnnotations
        String sampleId="sample"
    }

    if(defined(assemblyGfa)){
        call utils.extractAmbFromGFA {
            input:
            gfa=assemblyGfa
        }
    }
    if(!defined(assemblyGfa) && defined(assemblyPhase1) && defined(assemblyPhase0)){
        call utils.extractAmbFromGFAse {
            input:
            phase0=assemblyPhase0,
            phase1=assemblyPhase1,
            unphased=assemblyUnphased
        }
    }

    if(!defined(callVcf)){
        File asm0 = select_first([extractAmbFromGFA.asm0, extractAmbFromGFAse.asm0])
        File asm1 = select_first([extractAmbFromGFA.asm1, extractAmbFromGFAse.asm1])
        
        call dipdiff_t.dipdiff_t {
            input:
            ctgsPat=asm0,
            ctgsMat=asm1,
            reference=reference,
            vntrAnnotations=select_first([vntrAnnotations])
        }
    }

    File vcf_to_eval = select_first([callVcf, dipdiff_t.dipdiffUnphasedVcf])
    
    if(defined(truthVcf) && defined(confidentRegions)){
        call indexVcfFasta {
            input:
            call_vcf=vcf_to_eval,
            truth_vcf=truthVcf,
            reference_fa=reference
        }
        
        call truvariEvaluate {
            input:
            call_vcf=vcf_to_eval,
            call_vcf_index=indexVcfFasta.call_vcf_index,
            truth_vcf=truthVcf,
            truth_vcf_index=indexVcfFasta.truth_vcf_index,
            confident_bed=confidentRegions,
            reference_fa=reference,
            reference_fa_index=indexVcfFasta.reference_index
        }
    }
    
    output {
		File? vcfUnphased = dipdiff_t.dipdiffUnphasedVcf
		File? vcfPhased = dipdiff_t.dipdiffPhasedVcf
        File? evalSummary = truvariEvaluate.summary
        File? evalTarball = truvariEvaluate.output_tarball
    }
}

##
#### TASKS
##

task indexVcfFasta {
    input {
        File call_vcf
        File? truth_vcf
        File reference_fa
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="quay.io/biocontainers/samtools:1.16.1--h6899075_1"
        Int disk_size = 5 * round(size(call_vcf, 'G') + size(reference_fa, 'G')) + 20
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

        ln -s ~{call_vcf} call.vcf.gz
        ln -s ~{truth_vcf} truth.vcf.gz
        ln -s ~{reference_fa} ref.fa
        
        samtools faidx ref.fa
        tabix -p vcf call.vcf.gz
        tabix -p vcf truth.vcf.gz
    >>>
    output {
        File reference_index = "ref.fa.fai"
        File call_vcf_index = "call.vcf.gz.tbi"
        File truth_vcf_index = "truth.vcf.gz.tbi"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}

task truvariEvaluate {
    input {
        File call_vcf
        File call_vcf_index
        File? truth_vcf
        File? truth_vcf_index
        File? confident_bed
        File reference_fa
        File reference_fa_index
        String other_args = "-r 2000"
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="quay.io/biocontainers/truvari:3.5.0--pyhdfd78af_0"
        Int disk_size = 5 * round(size(call_vcf, 'G') + size(reference_fa, 'G')) + 20
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

        ln -s ~{call_vcf} call.vcf.gz
        ln -s ~{truth_vcf} truth.vcf.gz
        ln -s ~{call_vcf_index} call.vcf.gz.tbi
        ln -s ~{truth_vcf_index} truth.vcf.gz.tbi
        ln -s ~{reference_fa} ref.fa
        ln -s ~{reference_fa_index} ref.fa.fai

        truvari bench --no-ref a -b truth.vcf.gz -c call.vcf.gz -f ref.fa ~{other_args} --includebed ~{confident_bed} -o truvari_output/
        cp truvari_output/summary.txt summary.json
        tar -czvf truvari_output.tar.gz truvari_output
    >>>
    output {
        File summary = "summary.json"
        File output_tarball = "truvari_output.tar.gz"
    }
    runtime {
        docker: dockerImage
        cpu: threadCount
        disks: "local-disk " + disk_size + " SSD"
        memory: memoryGB + "GB"
    }
}


version 1.0

import "../tasks/gfase_utils.wdl" as utils
import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/master/QC/wdl/tasks/merqury.wdl" as merqury_t
import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/master/QC/wdl/tasks/meryl.wdl" as meryl_t

workflow QVtrio {

    input {
        File? assemblyGfa
        File? assemblyPhase0
        File? assemblyPhase1
        File? assemblyUnphased
        Array[File]? sampleReadsILM
        Array[File]? sampleReadsILMcram
        Array[File]? paternalReadsILM
        Array[File]? paternalReadsILMcram
        Array[File]? maternalReadsILM
        Array[File]? maternalReadsILMcram
        File? referenceFasta
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
    
    File asm0 = select_first([extractAmbFromGFA.asm0, extractAmbFromGFAse.asm0])
    File asm1 = select_first([extractAmbFromGFA.asm1, extractAmbFromGFAse.asm1])

    ## sample reads
    if(!defined(sampleReadsILM) && defined(sampleReadsILMcram) && defined(referenceFasta)){
        scatter (cram_file in select_first([sampleReadsILMcram])){
            call utils.convertCRAMtoFASTQ {
                input:
                cram_file=cram_file,
                reference_fa=select_first([referenceFasta])
            }
        }
    }

    Array[File] readFiles = select_first([sampleReadsILM, convertCRAMtoFASTQ.fastq_file])

    ## paternal reads
    if(!defined(paternalReadsILM) && defined(paternalReadsILMcram) && defined(referenceFasta)){
        scatter (cram_file in select_first([paternalReadsILMcram])){
            call utils.convertCRAMtoFASTQ as convertCRAMtoFASTQpaternal {
                input:
                cram_file=cram_file,
                reference_fa=select_first([referenceFasta])
            }
        }
    }

    Array[File] paternalFiles = select_first([paternalReadsILM, convertCRAMtoFASTQpaternal.fastq_file])

    ## maternal reads
    if(!defined(maternalReadsILM) && defined(maternalReadsILMcram) && defined(referenceFasta)){
        scatter (cram_file in select_first([maternalReadsILMcram])){
            call utils.convertCRAMtoFASTQ as convertCRAMtoFASTQmaternal {
                input:
                cram_file=cram_file,
                reference_fa=select_first([referenceFasta])
            }
        }
    }

    Array[File] maternalFiles = select_first([maternalReadsILM, convertCRAMtoFASTQmaternal.fastq_file])

    ## yak
    call utils.yakCount {
        input:
        readFiles=readFiles,
        sampleName=sampleId
    }

    call utils.yakCount as yakCountMaternal {
        input:
        readFiles=maternalFiles,
        sampleName=select_first(prefix("mat_", [sampleId]))
    }

    call utils.yakCount as yakCountPaternal {
        input:
        readFiles=paternalFiles,
        sampleName=select_first(prefix("pat_", [sampleId]))
    }

    call yakSwitchErr {
        input:
        yakCountFile=yakCount.yakCount,
        yakCountFileMaternal=yakCountMaternal.yakCount,
        yakCountFilePaternal=yakCountPaternal.yakCount,
        assemblyFastaHap1=asm0,
        assemblyFastaHap2=asm1,
        sampleName=sampleId
    }

    ## merqury
    call meryl_t.runMeryl as meryl {
        input:
        sampleReadsILM=readFiles,
        maternalReadsILM=maternalFiles,
        paternalReadsILM=paternalFiles,
        referenceFasta=referenceFasta
    }
    call merqury_t.merqury as merqury {
        input:
        assemblyFasta=asm0,
        altHapFasta=asm1,
        kmerTarball=meryl.sampleMerylDB,
        matKmerTarball=meryl.maternalHapmer,
        patKmerTarball=meryl.paternalHapmer
    }
    call extractMerquryFiles {
        input:
        merquryTarball=merqury.outputTarball
    }
    
    output {
        File yakSwitchErrSummary = yakSwitchErr.outputSummary
        File yakSwitchErrTarball = yakSwitchErr.outputTarball
        File merquryTarball = merqury.outputTarball
        File merquryQv = extractMerquryFiles.qv
        File merquryCompleteness = extractMerquryFiles.completeness
        File merqurySwitches = extractMerquryFiles.switches
    }
}

##
## TASKS
##

task yakSwitchErr {
    input{
        File yakCountFile
        File yakCountFileMaternal
        File yakCountFilePaternal
        File assemblyFastaHap2
        File assemblyFastaHap1
        String sampleName
        String genomeSize = "3.2g"
        String minSequenceLength = "100k"
        # runtime configurations
        Int memSizeGB=128
        Int threadCount=16
        Int diskSizeGB= 10 * round(size(yakCountFile, 'G') + size(yakCountFileMaternal, 'G') + size(yakCountFilePaternal, 'G') + size(assemblyFastaHap1, 'G') + size(assemblyFastaHap2, 'G')) + 50
        String dockerImage="juklucas/hpp_yak:latest"
    }
    command <<<
        set -eux -o xtrace

        # Computing error rates
        yak trioeval -t ~{threadCount} ~{yakCountFilePaternal} ~{yakCountFileMaternal} ~{assemblyFastaHap1} > ~{sampleName}.hap1.yak.switch-error.txt
        yak trioeval -t ~{threadCount} ~{yakCountFilePaternal} ~{yakCountFileMaternal} ~{assemblyFastaHap2} > ~{sampleName}.hap2.yak.switch-error.txt

        # condense
        SUMMARY=~{sampleName}.summary.txt
        echo "# hap1 switch" >>$SUMMARY
        tail -n3 ~{sampleName}.hap1.yak.switch-error.txt >>$SUMMARY
        echo "# hap2 switch" >>$SUMMARY
        tail -n3 ~{sampleName}.hap2.yak.switch-error.txt >>$SUMMARY

        # tar
        tar czvf ~{sampleName}.yak-switcherror.tar.gz *txt
    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File outputSummary = "~{sampleName}.summary.txt"
        File outputTarball = "~{sampleName}.yak-switcherror.tar.gz"
    }
}

task extractMerquryFiles {
    input{
        File merquryTarball
        String dockerImage="meredith705/gfase:latest"
    }
    command <<<
        set -eux -o xtrace

        tar -xzvf ~{merquryTarball}

        cat *.switches.txt > switches.txt
    >>>

    runtime {
        docker: dockerImage
        memory: "4 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
        preemptible: 1
    }

    output {
        File qv = glob("*.merqury.qv")[0]
        File completeness = glob("*.merqury.completeness.stats")[0]
        File switches = "switches.txt"
    }
}


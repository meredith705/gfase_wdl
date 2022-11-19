version 1.0

import "../tasks/gfase_utils.wdl" as utils

workflow QVnontrio {

    input {
        File? assemblyGfa
        File? assemblyPhase0
        File? assemblyPhase1
        File? assemblyUnphased
        Array[File] sampleReadsILM
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

    call yakQV {
        input:
        readFiles=sampleReadsILM,
        assemblyFastaHap1=asm0,
        assemblyFastaHap2=asm1,
        sampleName=sampleId
    }
            
    output {
        File yakSummary = yakQV.outputSummary
        File yakTarball = yakQV.outputTarball
    }
}

##
## TASKS
##

task yakQV {
    input{
        Array[File] readFiles
        File assemblyFastaHap2
        File assemblyFastaHap1
        String sampleName
        Int bloomSize=37
        String genomeSize = "3.2g"
        String minSequenceLength = "100k"
        # runtime configurations
        Int memSizeGB=128
        Int threadCount=16
        Int diskSizeGB= 10 * round(size(readFiles, 'G') + size(assemblyFastaHap1, 'G')) + 50
        String dockerImage="juklucas/hpp_yak:latest"
    }
    command <<<
        set -eux -o xtrace

        # Kmer counting with https://github.com/lh3/yak.
        yak count -t~{threadCount} -b~{bloomSize} -o ~{sampleName}.yak <(zcat ~{sep=" " readFiles}) <(zcat ~{sep=" " readFiles})

        # QV
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleName}.yak ~{assemblyFastaHap1} > ~{sampleName}.hap1.yak.qv.txt
        yak qv -t ~{threadCount} -p -K ~{genomeSize} -l ~{minSequenceLength} ~{sampleName}.yak ~{assemblyFastaHap2} > ~{sampleName}.hap2.yak.qv.txt

        # condense
        SUMMARY=~{sampleName}.summary.txt
        echo "# hap1 qv" >> $SUMMARY
        tail -n4 ~{sampleName}.hap1.yak.qv.txt >>$SUMMARY
        echo "# hap2 qv" >>$SUMMARY
        tail -n4 ~{sampleName}.hap2.yak.qv.txt >>$SUMMARY

        # tar
        tar czvf ~{sampleName}.yak-qc.tar.gz *txt
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
        File outputTarball = "~{sampleName}.yak-qc.tar.gz"
    }
}

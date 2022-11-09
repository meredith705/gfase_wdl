version 1.0

import "../tasks/whatshap.wdl" as whatshap_t
import "https://raw.githubusercontent.com/human-pangenomics/hpp_production_workflows/f8b09cb7b02729e57b5960665358e081ecfa6af6/QC/wdl/tasks/dipcall.wdl" as dipcall_t
import "../tasks/asmgene.wdl" as asmgene_t

workflow dipcallWhatshapAsmgene {

    input {
        File assemblyPat
        File assemblyMat
        File? assemblyUnphased
        File? truthVcf
        File reference
        String sampleId="sample"
        File asmGenesFasta
        File asmgeneReferenceFasta

    }

    if (defined(truthVcf)){
        call dipcall_t.dipcall as dipcall {
            input:
                assemblyFastaPat=assemblyPat,
                assemblyFastaMat=assemblyMat,
                referenceFasta=reference
        }

        call addPhaseSetToVCF {
            input:
                gzippedVcf=dipcall.outputVCF,
                outputIdentifier=sampleId,
        }

        call whatshap_t.whatshapAnalysis {
            input:
                queryVcf=addPhaseSetToVCF.phasettedVcf,
                truthVcf=truthVcf
        }

        call whatshap_t.coalesceResults {
            input:
                tarballs=[whatshapAnalysis.outputTarball],
                outputIdentifier=sampleId
        }
    }

    if (defined(assemblyUnphased)){
        call asmgene_t.asmgene as asmgenePatunp {
            input:
                assemblyFasta=assemblyPat,
                unphasedFasta=assemblyUnphased,
                genesFasta=asmGenesFasta,
                referenceFasta=asmgeneReferenceFasta
        }

        call asmgene_t.asmgene as asmgeneMatunp {
            input:
                assemblyFasta=assemblyMat,
                unphasedFasta=assemblyUnphased,
                genesFasta=asmGenesFasta,
                referenceFasta=asmgeneReferenceFasta
        }
    }

    if (!defined(assemblyUnphased)){
        call asmgene_t.asmgene as asmgenePat {
            input:
                assemblyFasta=assemblyPat,
                genesFasta=asmGenesFasta,
                referenceFasta=asmgeneReferenceFasta
        }

        call asmgene_t.asmgene as asmgeneMat {
            input:
            assemblyFasta=assemblyMat,
            genesFasta=asmGenesFasta,
            referenceFasta=asmgeneReferenceFasta
        }
    }

    File asmgeneGeneStatsPat = select_first([asmgenePat.geneStats,asmgenePatunp.geneStats] )
    File asmgeneGPStatsPat = select_first([asmgenePat.perGeneStats,asmgenePatunp.perGeneStats])
    File asmgeneGeneStatsMat = select_first([asmgeneMat.geneStats,asmgeneMatunp.geneStats ])
    File asmgeneGPStatsMat = select_first([asmgeneMat.perGeneStats,asmgeneMatunp.perGeneStats])


    output {
        File? whatshapTarball = coalesceResults.outputTarball
        File? whatshapReport = coalesceResults.fullOutput
        File asmgeneStatsPat = asmgeneGeneStatsPat
        File asmgenePGStatsPat = asmgeneGPStatsPat
        File asmgeneStatsMat = asmgeneGeneStatsMat
        File asmgenePGStatsMat = asmgeneGPStatsMat
		File? dipcallTarball = dipcall.outputTarball
		File? dipcallVCF = addPhaseSetToVCF.phasettedVcf
		File? dipcallBED = dipcall.outputBED
    }
}


task addPhaseSetToVCF {

    input {
        File gzippedVcf
        String outputIdentifier
        Int threadCount = 1
        Int memoryGB = 4
        String dockerImage="tpesout/whatshap:latest"
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

        OUTPUT=$(basename ~{gzippedVcf} | sed 's/.vcf.gz$//').PS.vcf

        zcat ~{gzippedVcf} | grep "^##" >$OUTPUT
        echo '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase Set Identifier">' >>$OUTPUT
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{outputIdentifier}" >>$OUTPUT

        zcat ~{gzippedVcf} | grep -v "^#" | sed 's/AD/AD:PS/' | sed 's/$/:1/' >>$OUTPUT

        bgzip $OUTPUT

    >>>

    output {
        File phasettedVcf = glob("*.PS.vcf.gz")[0]
    }

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: memoryGB + "GB"
    }
}


task makeSummaryStatFile {
    input {
        File dipWhatshapFull
        File whatsHapSWbed
        String assemblyID
        Int memSizeGB = 128
        Int threadCount = 1
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
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        # combine output files into one summary file
        # write to outputfile
        SUMMARY_FILE=`basename ~{assemblyID}`.summary.txt
        echo "dipcall/whatshap all_switch_rate:" >> $SUMMARY_FILE

        awk '{print $2,"\t",$32,"\t",$33,"\t",$37}' ~{dipWhatshapFull} >> $SUMMARY_FILE
    >>>

    runtime {
        docker: dockerImage
        cpu: threadCount
        memory: memSizeGB + "GB"
    }
}


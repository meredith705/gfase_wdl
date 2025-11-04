version 1.0

workflow runAsmgene {
	input{
        File assemblyFa
        File genesFa
        String sample
        File? unphasedFa
        File? referenceFa
        File? genesToReferencePaf
        Float asmgene_identity = 0.97
        # runtime configurations
        Int threadCount = 32
        Int memSizeGB = 64
        Int diskSizeGB = 32
        String dockerImage = "tpesout/hpp_minimap2:latest"
    }

    if (defined(unphasedFa)){
        call asmgene as asmgeneAsmUnp {
            input:
                sample=sample,
                assemblyFasta=assemblyFa,
                unphasedFasta=unphasedFa,
                genesFasta=genesFa,
                genesToReferencePaf=genesToReferencePaf,
                referenceFasta=referenceFa
        }

    }

    if (!defined(unphasedFa)){
        call asmgene as asmgeneAsm {
            input:
                sample=sample,
                assemblyFasta=assemblyFa,
                genesFasta=genesFa,
                genesToReferencePaf=genesToReferencePaf,
                referenceFasta=referenceFa
        }

    }

    File asmgeneGeneStatsAsm = select_first([asmgeneAsm.geneStats,asmgeneAsmUnp.geneStats] )
    File asmgeneGPStatsAsm = select_first([asmgeneAsm.perGeneStats,asmgeneAsmUnp.perGeneStats])


    output {
        File geneStats = asmgeneGeneStatsAsm
        File perGeneStats = asmgeneGPStatsAsm
        File? refGenePaf = asmgeneAsm.refPaf
        File? refGenePafup = asmgeneAsmUnp.refPaf
        File? asmGenePaf = asmgeneAsm.asmPaf
        File? asmGenePafup = asmgeneAsmUnp.asmPaf
    }


}

task asmgene {
    input{
        File assemblyFasta
        File genesFasta
        String sample
        String ref_name = "CHM13"
        File? unphasedFasta
        File? referenceFasta
        File? genesToReferencePaf
        Float asmgene_identity = 0.97
        Boolean outputRefPaf = false
        Boolean outputAsmPaf = false
        # runtime configurations
        Int threadCount = 32
        Int memSizeGB = 64
        Int diskSizeGB = 32
        String dockerImage = "tpesout/hpp_minimap2:latest"
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

        # name prep
        FILENAME=$(basename -- "~{assemblyFasta}" | sed 's/.gz$//' )
        PREFIX="${FILENAME%.*}"

        # determine if we need to align genes to reference
        if [[ -f "~{genesToReferencePaf}" ]] ; then
            ln -s ~{genesToReferencePaf} genesToRef.paf
        elif [[ -f "~{referenceFasta}" ]] ; then
            minimap2 -cx splice:hq -t ~{threadCount} ~{referenceFasta} ~{genesFasta} > genesToRef.paf
        else
            echo "Either referenceFasta or genesToReferencePaf needs to be supplied"
            exit 1
        fi

        cp ~{assemblyFasta} asm.fa
        if [[ -f "~{unphasedFasta}" ]]; then
            cat ~{unphasedFasta} >> asm.fa
            PREFIX=$PREFIX.unp
        fi

        # Aligning genes to assembly
        minimap2 -cx splice:hq -t ~{threadCount} asm.fa ~{genesFasta} > $PREFIX.paf

        # Computing statistics for gene completeness
        paftools.js asmgene -i~{asmgene_identity} -a genesToRef.paf $PREFIX.paf > $PREFIX.~{asmgene_identity}.gene_stats.txt

        # Output the stats for each gene to a file 
        paftools.js asmgene -i~{asmgene_identity} -e -a genesToRef.paf $PREFIX.paf > $PREFIX.~{asmgene_identity}.per_gene_stats.txt

        # determine if ref alignment is output
        if [[ "~{outputRefPaf}" == "true" ]]
        then
            # ref name prep
            REFFILENAME=$(basename -- "~{referenceFasta}" | sed 's/.gz$//' )
            REFPREFIX="${REFFILENAME%.*}"

            mv genesToRef.paf "~{ref_name}_genes.paf"
        fi

        # determine if asm alignment is output
        if [[ "~{outputAsmPaf}" == "true" ]]
        then
            mv "${PREFIX}.paf" "~{sample}_asmgenes.paf"
        fi

    >>>

    runtime {
        docker: dockerImage
        memory: memSizeGB + " GB"
        cpu: threadCount
        disks: "local-disk " + diskSizeGB + " SSD"
        preemptible: 1
    }

    output {
        File geneStats = glob("*.gene_stats.txt")[0]
        File perGeneStats = glob("*.per_gene_stats.txt")[0]
        File? refPaf = ref_name + "_genes.paf"
        File? asmPaf = sample + "_asmgenes.paf" 
    }
}


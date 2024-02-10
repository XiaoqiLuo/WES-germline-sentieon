version 1.0

workflow Annotate {
    input {
        File vcfFile
        File annovar
        Int NUM_THREAD = 20
        String MEMORY
        String DISK
    }
    output {
        File avinput = FormatAndAnnotate.avinput
        File out_csv = FormatAndAnnotate.out_csv
        
    }
    call FormatAndAnnotate {
        input:
            annovar = annovar,
            vcfFile = vcfFile,
            NUM_THREAD = NUM_THREAD,
            MEMORY = MEMORY,
            DISK = DISK
    }
}

task FormatAndAnnotate {
    input {
        File annovar
        File vcfFile
        Int NUM_THREAD = 20
        String MEMORY
        String DISK
    }

    command {
        set -e pipefail
        cp ${annovar} .
        tar -xzvf annovar.tar.gz
        CURRENT_DIR=$PWD
        cd ./annovar
        perl convert2annovar.pl -format vcf4old ${vcfFile} > $CURRENT_DIR/temp.avinput
        perl table_annovar.pl $CURRENT_DIR/temp.avinput humandb/ -buildver hg38 -out $CURRENT_DIR/output -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring . -csvout -polish -xref example/gene_xref.txt
    }


    output {
        File avinput = "temp.avinput"
        File out_csv = "output.hg38_multianno.csv"
    }

    runtime {
        docker: "registry-vpc.miracle.ac.cn/gznl/wes-annotate:v1.0.0"
        continueOnReturnCode: 0
        cpu: "${NUM_THREAD}"
        memory:"${MEMORY}" 
        disk: "${DISK}"
    }
}

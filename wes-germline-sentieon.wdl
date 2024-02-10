version 1.0

workflow DNAseq {
    
    # 流程输入文件和参数
    input {
        File fastq1
        File fastq2
        String sample_name
        String read_group_info
        File ref_fasta
        File interval
        File model
        File adapter
        Array[File] ref_fasta_indexes
        File dbsnp_vcf
        File dbsnp_index
        Array[File] bqsr_vcfs
        Array[File] bqsr_vcf_indexes
        String gotc_docker
        String sentieon_license_server
    }
    
    # 调用子任务
    call SentieonFastqToVcf {
        input:
            fastq1=fastq1,
            fastq2=fastq2,
            ref_fasta=ref_fasta,
            ref_fasta_indexes=ref_fasta_indexes,
            interval=interval,
            dbsnp_vcf=dbsnp_vcf,
            model=model,
            adapter=adapter,
            dbsnp_index=dbsnp_index,
            bqsr_vcfs=bqsr_vcfs,
            bqsr_vcf_indexes=bqsr_vcf_indexes,
            sample_name=sample_name,
            read_group_info=read_group_info,
            docker_image = gotc_docker,
            sentieon_license_server=sentieon_license_server
    }
    
    # 流程分析结果
    output {
        File recaled_bam = SentieonFastqToVcf.recaled_bam
        File recaled_bam_index = SentieonFastqToVcf.recaled_bam_index
        File bqsr_table = SentieonFastqToVcf.bqsr_table
        File gvcf = SentieonFastqToVcf.gvcf
        File gvcf_index = SentieonFastqToVcf.gvcf_index
        File aln_metrics = SentieonFastqToVcf.aln_metrics
        File dedup_metrics = SentieonFastqToVcf.dedup_metrics
        File gc_metrics = SentieonFastqToVcf.gc_metrics
        File gc_summary_metrics = SentieonFastqToVcf.gc_summary_metrics
        File is_metrics = SentieonFastqToVcf.is_metrics
        File qd_metrics = SentieonFastqToVcf.qd_metrics
    }
}

task SentieonFastqToVcf {
    input {
        File fastq1
        File fastq2
        File ref_fasta
        Array[File] ref_fasta_indexes
        File dbsnp_vcf
        File adapter
        File dbsnp_index
        File interval
        File model
        Array[File] bqsr_vcfs
        Array[File] bqsr_vcf_indexes
        String sample_name
        String read_group_info
        ## Runtime parameters
        String sentieon_license_server
        String docker_image
        Int NUM_THREAD = 20
        String MEMORY = "64 GB"
        String DISK = "500 GB"
        ## Extra algo parameters
        Int call_conf = 30
        String genotype_model = "Optional parameters: coalescent, multinomial. default multinomial"
    }

    Map[String,String] genotype_model_dict = {"Optional parameters: coalescent, multinomial. default multinomial": "multinomial",
                                              "multinomial":"multinomial",
                                              "coalescent": "coalescent"}
    String genotype_model_string = genotype_model_dict[genotype_model]

    String bqsr_table = sample_name + ".recal_data.table"
    String out_gvcf = sample_name + ".gvcf.gz"
    String out_gvcf_idx = sample_name + ".gvcf.gz.tbi"

    String dollar = "$"
    
    # 工具运行命令
    command <<<
        
        export SENTIEON_LICENSE=~{sentieon_license_server}

        ## Step0
        fastp -i ~{fastq1} -I ~{fastq2} -o ~{sample_name}_trim_1.fq.gz -O ~{sample_name}_trim_2.fq.gz -g -q 5 -u 50 -n 15 -w 1 \
        -j ~{sample_name}.json -h ~{sample_name}.html

        clean_fq1=~{sample_name}_trim_1.fq.gz
        clean_fq2=~{sample_name}_trim_2.fq.gz
     
        ## Step1
        ( sentieon bwa mem -R ~{read_group_info} -K 100000000 -t $(nproc) ~{ref_fasta} $clean_fq1 $clean_fq2 || echo -n 'error' ) \
            | sentieon util sort -i - -r ~{ref_fasta} -t $(nproc) -o ~{sample_name}.sorted.bam --sam2bam

        next_bam=~{sample_name}.sorted.bam

        # MeanQualityByCycle
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{sample_name}.sorted.bam  \
         --algo MeanQualityByCycle ~{sample_name}.mq_metrics.txt \
         --algo QualDistribution ~{sample_name}.qd_metrics.txt \
         --algo GCBias --summary ~{sample_name}.gc_summary_metrics.txt ~{sample_name}.gc_metrics.txt \
         --algo AlignmentStat ~{sample_name}.aln_metrics.txt \
         --algo InsertSizeMetricAlgo ~{sample_name}.is_metrics.txt

        # Step2
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i $next_bam \
         --algo LocusCollector \
         ~{sample_name}.score.txt.gz

        # Step3
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i $next_bam \
            --algo Dedup \
            --score_info ~{sample_name}.score.txt.gz \
            --metrics ~{sample_name}.dedup_metrics.txt \
            --rmdup \
            ~{sample_name}.dedup.bam
        next_bam=~{sample_name}.dedup.bam
        
        # Step4
        touch ~{bqsr_table}
        
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i $next_bam  \
           --algo QualCal \
           ~{"-k "+ dbsnp_vcf} \
           -k ~{sep=" -k " bqsr_vcfs} \
           ~{bqsr_table}

        bqsr_string="--read_filter QualCalFilter,table=~{bqsr_table},prior=-1.0,indel=false,levels=10/20/30,min_qual=6"
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i $next_bam $bqsr_string \
            --algo ReadWriter ~{sample_name}.recaled.bam
        next_bam=~{sample_name}.recaled.bam

        #Step5
        sentieon driver -r ~{ref_fasta} -t $(nproc) --interval ~{interval} \
            --interval_padding 100 -i $next_bam -q ~{bqsr_table} --algo Haplotyper --emit_mode gvcf -d ~{dbsnp_vcf} ~{out_gvcf}

    >>>

    runtime {
        docker: docker_image
        cpu: "${NUM_THREAD}"
        memory:"${MEMORY}" 
        disk: "${DISK}"
    }
    
    output {
        File recaled_bam = "~{sample_name}.recaled.bam"
        File recaled_bam_index = "~{sample_name}.recaled.bam.bai"
        File bqsr_table = bqsr_table
        File gvcf = "~{out_gvcf}"
        File gvcf_index = "~{out_gvcf}.tbi"
        File aln_metrics = "~{sample_name}.aln_metrics.txt"
        File dedup_metrics = "~{sample_name}.dedup_metrics.txt"
        File gc_metrics = "~{sample_name}.gc_metrics.txt"
        File gc_summary_metrics = "~{sample_name}.gc_summary_metrics.txt"
        File is_metrics = "~{sample_name}.is_metrics.txt"
        File qd_metrics = "~{sample_name}.qd_metrics.txt"
    }
}


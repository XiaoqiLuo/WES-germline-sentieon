version 1.0

workflow DNAseq {
    
    # 流程输入文件和参数
    input {
        File tumor_fastq1
        File tumor_fastq2
        File normal_fastq1
        File normal_fastq2
        String TUMOR_SM
        String NORMAL_SM
        String read_group_info_NORMAL
        String read_group_info_TUMOR
        String sample_name
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
            tumor_fastq1=tumor_fastq1,
            tumor_fastq2=tumor_fastq2,
            normal_fastq1=normal_fastq1,
            normal_fastq2=normal_fastq2,
            TUMOR_SM=TUMOR_SM,
            NORMAL_SM=NORMAL_SM,
            read_group_info_NORMAL=read_group_info_NORMAL,
            read_group_info_TUMOR=read_group_info_TUMOR,
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
            # read_group_info=read_group_info,
            docker_image = gotc_docker,
            sentieon_license_server=sentieon_license_server
    }
    
    # 流程分析结果
  output {
    File tumor_sorted_bam = SentieonFastqToVcf.tumor_sorted_bam
    File normal_sorted_bam = SentieonFastqToVcf.normal_sorted_bam
    File normal_mq_metrics = SentieonFastqToVcf.normal_mq_metrics
    File normal_qd_metrics = SentieonFastqToVcf.normal_qd_metrics
    File normal_gc_metrics = SentieonFastqToVcf.normal_gc_metrics
    File normal_gc_summary_metrics = SentieonFastqToVcf.normal_gc_summary_metrics
    File normal_aln_metrics = SentieonFastqToVcf.normal_aln_metrics
    File normal_is_metrics = SentieonFastqToVcf.normal_is_metrics
    File tumor_dedup_bam = SentieonFastqToVcf.tumor_dedup_bam
    File tumor_dedup_metrics = SentieonFastqToVcf.tumor_dedup_metrics
    File normal_dedup_bam = SentieonFastqToVcf.normal_dedup_bam
    # File normal_dedup_metrics = SentieonFastqToVcf.normal_dedup_metrics
    # File normal_coverage_metrics = SentieonFastqToVcf.normal_coverage_metrics
    # File tumor_coverage_metrics = SentieonFastqToVcf.tumor_coverage_metrics
    File tumor_recal_data_table = SentieonFastqToVcf.tumor_recal_data_table
    File tumor_recal_data_table_post = SentieonFastqToVcf.tumor_recal_data_table_post
    File tumor_recal_csv = SentieonFastqToVcf.tumor_recal_csv
    File tumor_recal_plots_pdf = SentieonFastqToVcf.tumor_recal_plots_pdf
    File normal_recal_data_table = SentieonFastqToVcf.normal_recal_data_table
    File normal_recal_data_table_post = SentieonFastqToVcf.normal_recal_data_table_post
    File normal_recal_csv = SentieonFastqToVcf.normal_recal_csv
    File normal_recal_plots_pdf = SentieonFastqToVcf.normal_recal_plots_pdf
    File tmp_vcf = SentieonFastqToVcf.tmp_vcf
    File orientation_bias = SentieonFastqToVcf.orientation_bias
    File final_vcf = SentieonFastqToVcf.final_vcf
    }

}

task SentieonFastqToVcf {
    input {
        File tumor_fastq1
        File tumor_fastq2
        File normal_fastq1
        File normal_fastq2
        File ref_fasta
        String NORMAL_SM
        String TUMOR_SM
        Array[File] ref_fasta_indexes
        File dbsnp_vcf
        File adapter
        File dbsnp_index
        File interval
        File model
        Array[File] bqsr_vcfs
        Array[File] bqsr_vcf_indexes
        String sample_name
        String read_group_info_NORMAL
        String read_group_info_TUMOR
        ## Runtime parameters
        String sentieon_license_server
        String docker_image
        Int NUM_THREAD = 30
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

        # /root/result-data/lxq/sentieon/sentieon-scripts-master/example_pipelines/somatic/TNseq/tumor_normal.sh
        ## Step1 tumor mapping
        ( sentieon bwa mem -R ~{read_group_info_TUMOR} -K 100000000 -Y -t $(nproc) ~{ref_fasta} ~{tumor_fastq1} ~{tumor_fastq2} || echo -n 'error' ) \
            | sentieon util sort -i - -r ~{ref_fasta} -t $(nproc) -o ~{TUMOR_SM}_tumor.sorted.bam --sam2bam

        ( sentieon bwa mem -R ~{read_group_info_NORMAL} -K 100000000 -Y -t $(nproc) ~{ref_fasta} ~{normal_fastq1} ~{normal_fastq2} || echo -n 'error' ) \
            | sentieon util sort -i - -r ~{ref_fasta} -t $(nproc) -o ~{NORMAL_SM}_normal.sorted.bam --sam2bam


        # Step2 Metrics for normal sample normal sample
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.sorted.bam  \
         --algo MeanQualityByCycle ~{NORMAL_SM}_normal.mq_metrics.txt \
         --algo QualDistribution ~{NORMAL_SM}_normal.qd_metrics.txt \
         --algo GCBias --summary ~{NORMAL_SM}_normal.gc_summary_metrics.txt ~{NORMAL_SM}_normal.gc_metrics.txt \
         --algo AlignmentStat ~{NORMAL_SM}_normal.aln_metrics.txt \
         --algo InsertSizeMetricAlgo ~{NORMAL_SM}_normal.is_metrics.txt

        # Step3a. Remove Duplicate Reads for tumor sample
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{TUMOR_SM}_tumor.sorted.bam \
         --algo LocusCollector \
         ~{TUMOR_SM}_tumor.score.txt.gz
        
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{TUMOR_SM}_tumor.sorted.bam \
            --algo Dedup \
            --score_info ~{TUMOR_SM}_tumor.score.txt.gz \
            --metrics ~{TUMOR_SM}_tumor.dedup_metrics.txt \
            --rmdup ~{TUMOR_SM}_tumor.dedup.bam
        
        cp ~{TUMOR_SM}_tumor.dedup.bam  ~{TUMOR_SM}.bam
        cp ~{NORMAL_SM}_normal.dedup.bam  ~{NORMAL_SM}.bam
        # Step3b. Remove Duplicate Reads for normal sample
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.sorted.bam \
            --algo LocusCollector \
            ~{NORMAL_SM}_normal.score.txt.gz
        
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.sorted.bam \
            --algo Dedup \
            --score_info ~{NORMAL_SM}_normal.score.txt.gz \
            --metrics ~{NORMAL_SM}_normal.dedup_metrics.txt \
            --rmdup ~{NORMAL_SM}_normal.dedup.bam

        # ******************************************
        # 4b. Coverage metrics for sample
        # ******************************************
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.dedup.bam \
            --algo CoverageMetrics normal_coverage_metrics 

        sentieon driver -r ~{ref_fasta} -t $(nproc) ~{TUMOR_SM}_tumor.dedup.bam \
            --algo CoverageMetrics tumor_coverage_metrics 


        # ******************************************
        # 5a. Base recalibration for tumor sample
        # ******************************************
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{TUMOR_SM}_tumor.dedup.bam \
           --algo QualCal \
           ~{"-k "+ dbsnp_vcf} \
           -k ~{sep=" -k " bqsr_vcfs} \
           ~{TUMOR_SM}_tumor_recal_data.table

        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{TUMOR_SM}_tumor.dedup.bam \
        -q ~{TUMOR_SM}_tumor_recal_data.table --algo QualCal \
           ~{"-k "+ dbsnp_vcf} \
           -k ~{sep=" -k " bqsr_vcfs}  ~{TUMOR_SM}_tumor_recal_data.table.post

        sentieon  driver -t $(nproc) --algo QualCal --plot --before \
            ~{TUMOR_SM}_tumor_recal_data.table --after ~{TUMOR_SM}_tumor_recal_data.table.post \
            ~{TUMOR_SM}_tumor_recal.csv

        sentieon plot QualCal -o ~{TUMOR_SM}_tumor_recal_plots.pdf ~{TUMOR_SM}_tumor_recal.csv

        # ******************************************
        # 5b. Base recalibration for normal sample
        # ******************************************
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.dedup.bam \
           --algo QualCal \
           ~{"-k "+ dbsnp_vcf} \
           -k ~{sep=" -k " bqsr_vcfs} \
           ~{NORMAL_SM}_normal_recal_data.table

        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.dedup.bam \
        -q ~{NORMAL_SM}_normal_recal_data.table --algo QualCal \
           ~{"-k "+ dbsnp_vcf} \
           -k ~{sep=" -k " bqsr_vcfs}  ~{NORMAL_SM}_normal_recal_data.table.post

        sentieon  driver -t $(nproc) --algo QualCal --plot --before \
            ~{NORMAL_SM}_normal_recal_data.table --after ~{NORMAL_SM}_normal_recal_data.table.post \
            ~{NORMAL_SM}_normal_recal.csv

        sentieon plot QualCal -o ~{NORMAL_SM}_normal_recal_plots.pdf ~{NORMAL_SM}_normal_recal.csv

        # ******************************************
        # 6. Somatic Variant Calling - TNhaplotyper2
        # ******************************************
        sentieon driver -r ~{ref_fasta} -t $(nproc) -i ~{NORMAL_SM}_normal.dedup.bam \
            -i ~{TUMOR_SM}_tumor.dedup.bam -q ~{NORMAL_SM}_normal_recal_data.table -q ~{TUMOR_SM}_tumor_recal_data.table \
                --algo TNhaplotyper2 \
                --tumor_sample ~{TUMOR_SM} \
                --normal_sample ~{NORMAL_SM} \
                output-tnhap2-tmp.vcf.gz \
                --algo OrientationBias --tumor_sample ~{TUMOR_SM} output-orientation 

        sentieon driver -r ~{ref_fasta} --algo TNfilter \
            -v output-tnhap2-tmp.vcf.gz --tumor_sample ~{TUMOR_SM} --normal_sample ~{NORMAL_SM} \
            --orientation_priors output-orientation output-tnhap2.vcf.gz 
    >>>

    runtime {
        docker: docker_image
        cpu: "${NUM_THREAD}"
        memory:"${MEMORY}" 
        disk: "${DISK}"
    }
    
    output {
    File tumor_sorted_bam = "~{TUMOR_SM}_tumor.sorted.bam"
    File normal_sorted_bam = "~{NORMAL_SM}_normal.sorted.bam"
    File normal_mq_metrics = "~{NORMAL_SM}_normal.mq_metrics.txt"
    File normal_qd_metrics = "~{NORMAL_SM}_normal.qd_metrics.txt"
    File normal_gc_metrics = "~{NORMAL_SM}_normal.gc_metrics.txt"
    File normal_gc_summary_metrics = "~{NORMAL_SM}_normal.gc_summary_metrics.txt"
    File normal_aln_metrics = "~{NORMAL_SM}_normal.aln_metrics.txt"
    File normal_is_metrics = "~{NORMAL_SM}_normal.is_metrics.txt"
    File tumor_dedup_bam = "~{TUMOR_SM}_tumor.dedup.bam"
    File tumor_dedup_metrics = "~{TUMOR_SM}_tumor.dedup_metrics.txt"
    File normal_dedup_bam = "~{NORMAL_SM}_normal.dedup.bam"
    File normal_dedup_metrics = "~{NORMAL_SM}_normal.dedup_metrics.txt"
    # File normal_coverage_metrics = "normal_coverage_metrics"
    # File tumor_coverage_metrics = "tumor_coverage_metrics"
    File tumor_recal_data_table = "~{TUMOR_SM}_tumor_recal_data.table"
    File tumor_recal_data_table_post = "~{TUMOR_SM}_tumor_recal_data.table.post"
    File tumor_recal_csv = "~{TUMOR_SM}_tumor_recal.csv"
    File tumor_recal_plots_pdf = "~{TUMOR_SM}_tumor_recal_plots.pdf"
    File normal_recal_data_table = "~{NORMAL_SM}_normal_recal_data.table"
    File normal_recal_data_table_post = "~{NORMAL_SM}_normal_recal_data.table.post"
    File normal_recal_csv = "~{NORMAL_SM}_normal_recal.csv"
    File normal_recal_plots_pdf = "~{NORMAL_SM}_normal_recal_plots.pdf"
    File tmp_vcf = "output-tnhap2-tmp.vcf.gz"
    File orientation_bias = "output-orientation"
    File final_vcf = "output-tnhap2.vcf.gz"
}

}


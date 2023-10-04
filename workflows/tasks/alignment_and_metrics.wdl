version 1.0

# Map the amplicon reads to the reference genome
task alignment_metrics {

    input {
        File genome_reference
        File genome_index
        File amplicons_fastq_gz
        String file_label
        String docker
        Int alignment_thread = 4
        Int sort_thread = 4
    }
 
    String log_level = "DEBUG"
    Int memory_mb = ceil(size(amplicons_fastq_gz, "GiB")) + 32
    Int disk_size_gb = ceil(size(amplicons_fastq_gz, "GiB")) + 15

    command <<<
        
        # initial quality check
        fastqc ~{amplicons_fastq_gz} -o .
        rm *_fastqc.html

        # align to reference genome and sort the output
        pbmm2 align \
        --log-level ~{log_level} \
        --log-file ~{file_label}_aligned_amplicon_cluster.log \
        --sort -j ~{alignment_thread} -J ~{sort_thread} \
        --preset HIFI \
        ~{genome_reference} ~{amplicons_fastq_gz} ~{file_label}_aligned_sorted.bam

        # Get coverage metrics
        samtools mpileup \
            --min-BQ 1 \
            -f ~{genome_reference} \
            -s ~{file_label}_aligned_sorted.bam > ~{file_label}_aligned_sorted.bam.mpileup

        samtools depth \
            -q 0 -Q 0 \
            ~{file_label}_aligned_sorted.bam > ~{file_label}_aligned_sorted.bam.depth
        
        samtools flagstat \
            ~{file_label}_aligned_sorted.bam >  ~{file_label}_aligned_final_sorted.flagstat

        samtools idxstat \
            ~{file_label}_aligned_sorted.bam >  ~{file_label}_aligned_final_sorted.idxstat
    >>>

    output {
        File fastqc_report = file_label + "_fastqc.zip"
        File aligned_sorted_bam = file_label + "_aligned_sorted.bam"
        File aligned_amplicon_cluster_log = file_label + "_aligned_amplicon_cluster.log"
        File bam_mpileup_report = file_label + "_aligned_sorted.bam.mpileup"
        File bam_depth_report = file_label + "_aligned_sorted.bam.depth"
        File flagstat_report = file_label + "_aligned_final_sorted.flagstat"
        File idxstat_report = file_label + "_aligned_final_sorted.idxstat"
    }

    parameter_meta {

    }

    runtime {
        docker: "~{docker}"
        memory: "~{memory_mb} MiB"
        disks: "local-disk ~{disk_size_gb} HDD"
        returnCodes: "*"
        continueOnReturnCode: true
    }
}
version 1.0

workflow main {

    String container_src = "pbaaenv2:latest"

    input {
        File reads_fastq_gz
        File guide_fasta
        File genome_ref
        String prefix
    }

    call pbmm_index{
        input: genome_reference = genome_ref, docker = container_src
    }
    call pbaa_cluster_var_call {
        input : pbaa_guide_fasta = guide_fasta, genome_reference = genome_ref, amplicons_fastq_gz = reads_fastq_gz, file_label = prefix, docker = container_src
    }
    call pbmm_alignment_metrics {
        input: genome_reference = genome_ref, genome_index = pbmm_index.index, amplicons_fastq_gz = reads_fastq_gz, file_label = prefix, docker = container_src
    }
    call vcfcons_variants_bam {
        input : genome_reference = genome_ref, genome_index = pbmm_index.index, bam_depth = pbmm_alignment_metrics.bam_depth_report, pbaa_vcf = pbaa_cluster_var_call.pbaa_vcf, file_label = prefix, docker = container_src
    }
}

task pbmm_index {
    input {
        File genome_reference
        String docker
    }

    command <<<
        pbmm2 index ~{genome_reference} ref.mmi --preset SUBREAD
    >>>

    output {
        File index = "ref.mmi"
    }

    runtime {
        docker: "~{docker}"
    }


}


# cluster using pbAA and variant calling
task pbaa_cluster_var_call {
    input {
        File pbaa_guide_fasta
        File genome_reference
        File amplicons_fastq_gz
        String file_label
        String docker
        Int min_cluster_read_count = 2
    }  

    String log_level = "DEBUG"

    command <<<

        set -euo pipefail
        
        fastq_gz_inp=~{file_label}.fastq.gz # replace ".fastq.gz" with ".hifi_reads.fastq.gz"
        
        fastq_inp=~{file_label}.fastq # .hifi_reads.fastq

        samtools faidx ~{pbaa_guide_fasta}

        samtools faidx ~{genome_reference}

        gunzip -c ~{amplicons_fastq_gz} > $fastq_inp

        samtools faidx --fastq $fastq_inp

        pbaa cluster \
            --log-level ~{log_level} \
            --log-file ~{file_label}_pbaa.log \
            --trim-ends 5 \
            --min-cluster-read-count ~{min_cluster_read_count} \
            ~{pbaa_guide_fasta} ${fastq_inp} ~{file_label}_pbaa

        # convert pbaa outcome to VCF
        python3 /scripts/consensusVariants.py \
            --runName ~{file_label}_pbaa \
            --prefix ~{file_label}_pbaa \
            --read_info ~{file_label}_pbaa_read_info.txt \
            --hifiSupport ${fastq_inp} \
            ~{genome_reference} \
            ~{file_label}_pbaa_passed_cluster_sequences.fasta > ~{file_label}_consensusVariants.log
        
        python3 /scripts/pbaa2vcf.py \
            --passOnly -s barcode \
            -o ~{file_label}.vcf \
            ~{file_label}_pbaa_alleles.csv \
            ~{file_label}_pbaa_variants.csv \
            ~{genome_reference} > ~{file_label}_pbaa2vcf.log
    >>>

    output {
        File passed_cluster_sequences = file_label + "_pbaa_passed_cluster_sequences.fasta"
        File failed_cluster_sequences = file_label + "_pbaa_failed_cluster_sequences.fasta"
        File read_info = file_label + "_pbaa_read_info.txt"
        File pbaa_run_log = file_label + "_pbaa.log"
        File pbaa_alleles = file_label + "_pbaa_alleles.csv"
        File pbaa_variants = file_label + "_pbaa_variants.csv"
        File pbaa_vcf = file_label + ".vcf"
    }

    parameter_meta {

    }

    runtime {
        docker: "~{docker}"
    }
}

# Map the amplicon reads to the reference genome
task pbmm_alignment_metrics {

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

    command <<<
        
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
    >>>

    output {
        File aligned_sorted_bam = file_label + "_aligned_sorted.bam"
        File aligned_amplicon_cluster_log = file_label + "_aligned_amplicon_cluster.log"
        File bam_mpileup_report = file_label + "_aligned_sorted.bam.mpileup"
        File bam_depth_report = file_label + "_aligned_sorted.bam.depth"
    }

    parameter_meta {

    }

    runtime {
        docker: "~{docker}"
    }
}


# get the final variants VCF
task vcfcons_variants_bam {
    input {
        File genome_reference
        File genome_index
        File bam_depth
        File pbaa_vcf
        String file_label
        String docker
        Int alignment_thread = 4
        Int sort_thread = 4
        Int min_coverage = 4
        Float min_alt_freq = 0.5
    }
    
    String log_level = "DEBUG"

    command <<<

        python3 /scripts/VCFCons.py \
            ~{genome_reference} ~{file_label} \
            --sample-name ~{file_label} \
            --min_coverage ~{min_coverage} \
            --min_alt_freq ~{min_alt_freq} \
            --vcf_type pbaa \
            --input_depth ~{bam_depth} \
            --input_vcf ~{pbaa_vcf} > ~{file_label}.vcfcons.log

        # Now map the consensus sequences to the reference genome
        pbmm2 align \
            --log-level ~{log_level} \
            --log-file ~{file_label}_consensus_reference_alignment.log \
            --sort -j ~{alignment_thread} -J ~{sort_thread} \
            --preset HIFI \
            ~{genome_reference} ~{file_label}.vcfcons.frag.fasta \
            ~{file_label}_vcfcons_aligned_sorted.bam
    >>>

    output {        
        File vcfcons_fasta = file_label + ".vcfcons.fasta"
        File vcfcons_frag_fasta = file_label + ".vcfcons.frag.fasta"
        File vcfcons_info = file_label + ".vcfcons.info.csv"
        File vcfcons_vcf = file_label + ".vcfcons.vcf" # VCF output from pbaa variant call
        File vcfcons_csv = file_label + ".vcfcons.variants.csv"
        File vcfcons_log = file_label + ".vcfcons.log"
        File consensus_reference_aligned_bam = file_label + "_vcfcons_aligned_sorted.bam"
        File consensus_reference_alignment_log = file_label + "_consensus_reference_alignment.log"
        # depth file of the input aligned to ref genome

    }

    parameter_meta {

    }

    runtime {
        docker: "~{docker}"
    }
}

# to do
# bam - alignment metrics: samtools flagstat, idxstat
# vcf - bcftools stats

# get the summary metrics
# task summary {

# }
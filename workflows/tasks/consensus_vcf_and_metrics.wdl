version 1.0

# get the final variants VCF using VCFCons
task consensus_variant_calling {
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
    # Int memory_mb = ceil(size(genome_index, "GiB")) + 32
    # Int disk_size_gb = ceil(size(genome_index, "GiB")) + 15

    command <<<

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{genome_index} genome_reference.fasta.mmi

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai
        
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
            genome_reference.fasta ~{file_label}.vcfcons.frag.fasta \
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
        memory: "32G"
        disks: "local-disk 25 HDD"
        returnCodes: "*"
        continueOnReturnCode: true
    }
}
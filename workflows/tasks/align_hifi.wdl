version 1.0

# align clustered hifi reads to reference and generate the final bam using pbaa bampaint
task HifiReadsAlign {
    
    input {
        File hifi_reads_fastq_gz
        File pbmm2_index
        String file_label
        Int subset
        String docker
    }  

    String log_level = "DEBUG"

    command <<<
        set -euo pipefail

        ln -s ~{pbmm2_index} genome_reference.mmi

        hifireads_file_base=$(basename ~{hifi_reads_fastq_gz} .fastq.gz)

        gunzip -c ~{hifi_reads_fastq_gz} > ${hifireads_file_base}.fastq

        if [[ ~{subset} -gt 0 && ~{subset} -lt 100 ]]
        then
            total_reads=$(awk '{s++}END{print s/4}' ${hifireads_file_base}.fastq)
            desired_percent_reads=$(awk -v total=$total_reads -v percent=~{subset} 'BEGIN{printf "%.0f", total*percent/100}')
            seqtk sample ${hifireads_file_base}.fastq ${desired_percent_reads} > ~{file_label}.hifi_reads.fastq
        else
            mv ${hifireads_file_base}.fastq ~{file_label}.hifi_reads.fastq
        fi
        
        samtools faidx --fastq ~{file_label}.hifi_reads.fastq

        pbmm2 align \
            --preset hifi \
            --sort \
            genome_reference.mmi \
            ~{file_label}.hifi_reads.fastq \
            ~{file_label}_raw_hifi_to_reference_alignment.bam \
            --rg "@RG\tID:~{file_label}\tSM:~{file_label}" \
            --log-level ~{log_level} \
            --log-file ~{file_label}_raw_hifi_to_reference_alignment.log
        
        seqkit stats -a -T ~{file_label}.hifi_reads.fastq > ~{file_label}_raw_hifi_reads_fastq_stats.tab
        
        samtools flagstat \
            ~{file_label}_raw_hifi_to_reference_alignment.bam > ~{file_label}_raw_hifi_to_reference_alignment_flagstat.txt
        
        samtools idxstat \
            ~{file_label}_raw_hifi_to_reference_alignment.bam > ~{file_label}_raw_hifi_to_reference_alignment_idx.txt

        printf "ref_contig\tref_contig_len\ttotal_mapped_reads\ttotal_unmapped_reads\n" | cat - ~{file_label}_raw_hifi_to_reference_alignment_idx.txt > ~{file_label}_raw_hifi_to_reference_alignment_idxstat.txt
    >>>

    output {
        File raw_hifi_to_reference_alignment_bam = file_label + "_raw_hifi_to_reference_alignment.bam"
        File raw_hifi_to_reference_alignment_index = file_label + "_raw_hifi_to_reference_alignment.bam.bai"
        File raw_hifi_to_reference_alignment_log = file_label + "_raw_hifi_to_reference_alignment.log"
        File raw_hifi_to_reference_alignment_flagstat = file_label + "_raw_hifi_to_reference_alignment_flagstat.txt"
        File raw_hifi_to_reference_alignment_idxstat = file_label + "_raw_hifi_to_reference_alignment_idxstat.txt"
        File raw_hifi_reads_fastq_stats = file_label + "_raw_hifi_reads_fastq_stats.tab"
    }

    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
version 1.0

# align consensus reads to reference
task alignConsensus {
    
    input {
        File genome_reference
        File pbmm2_index
        File passed_cluster_seq
        String file_label
        String docker
    }  

    String log_level = "DEBUG"

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{pbmm2_index} genome_reference.mmi
        
        pbmm2 align \
            --preset hifi \
            --sort \
            ~{genome_reference} \
            ~{passed_cluster_seq} \
            ~{file_label}_consensus_to_reference_alignment.bam \
            --rg "@RG\tID:~{file_label}\tSM:~{file_label}" \
            --log-level ~{log_level} \
            --log-file ~{file_label}_consensus_to_reference_alignment.log
        
        samtools flagstat \
            ~{file_label}_consensus_to_reference_alignment.bam > ~{file_label}_consensus_to_reference_alignment_flagstat.txt
        
        samtools idxstat \
            ~{file_label}_consensus_to_reference_alignment.bam > ~{file_label}_consensus_to_reference_alignment_idx.txt

        printf "ref_contig\tref_contig_len\ttotal_mapped_reads\ttotal_unmapped_reads\n" | cat - ~{file_label}_consensus_to_reference_alignment_idx.txt > ~{file_label}_consensus_to_reference_alignment_idxstat.txt
    >>>

    output {
        File consensus_to_reference_alignment_bam = file_label + "_consensus_to_reference_alignment.bam"
        File consensus_to_reference_alignment_bam_idx = file_label + "_consensus_to_reference_alignment.bam.bai"
        File consensus_to_reference_alignment_log = file_label + "_consensus_to_reference_alignment.log"
        File consensus_to_reference_alignment_flagstat = file_label + "_consensus_to_reference_alignment_flagstat.txt"
        File consensus_to_reference_alignment_idxstat = file_label + "_consensus_to_reference_alignment_idxstat.txt"
    }
    
    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
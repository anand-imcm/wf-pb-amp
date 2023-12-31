version 1.0

# call variants using the consensus to reference aligned bam
task variantCall {
    
    input {
        File consensus_to_ref_aligned_bam
        File consensus_to_ref_aligned_bam_index
        File genome_reference
        File clinvar
        File gff
        File bed
        String file_label
        String docker
    } 

    command <<<
        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{clinvar} clinvar.vcf.gz

        ln -s ~{consensus_to_ref_aligned_bam} ~{file_label}_consensus_to_ref_aligned.bam
        ln -s ~{consensus_to_ref_aligned_bam_index} ~{file_label}_consensus_to_ref_aligned.bam.bai

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai
        
        tabix -p vcf clinvar.vcf.gz
        
        python /scripts/variant_call.py \
            --inbam ~{file_label}_consensus_to_ref_aligned.bam \
            --ref genome_reference.fasta \
            --clinvar clinvar.vcf.gz \
            --gff ~{gff} \
            --target ~{bed} \
            --prefix ~{file_label}
    >>>

    output {
        File variant_summary = file_label + "_variant_summary.tsv"
        File variant_on_target_summary = file_label + "_variant_on_target_summary.tsv"
        Array[File] annotated_vcf = glob("*_annotated.vcf.gz")
        Array[File] raw_vcf = glob("*_raw.vcf.gz")
    }

    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
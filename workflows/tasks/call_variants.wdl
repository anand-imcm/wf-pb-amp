version 1.0

# align 
task varcall {
    
    input {
        File consensus_to_ref_aligned_bam
        File genome_reference
        File clinvar
        File gff
        File bed
        String file_label
    } 

    command <<<
        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{clinvar} clinvar.vcf.gz

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai
        
        tabix -p vcf clinvar.vcf.gz
        
        python /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/scripts/variant_call.py \
            --inbam ~{consensus_to_ref_aligned_bam} \
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
}
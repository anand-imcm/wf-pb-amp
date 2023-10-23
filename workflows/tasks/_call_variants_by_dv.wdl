version 1.0

# call variants using the consensus to reference aligned bam
task variantCallGDV {
    
    input {
        File consensus_to_ref_aligned_bam
        File consensus_to_ref_aligned_bam_index
        File genome_reference
        File clinvar
        File gff
        File bed
        String file_label
        String deepvariant_version = "1.5.0"
    } 

    command <<<
        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{clinvar} clinvar.vcf.gz

        ln -s ~{consensus_to_ref_aligned_bam} ~{file_label}_consensus_to_ref_aligned.bam
        ln -s ~{consensus_to_ref_aligned_bam_index} ~{file_label}_consensus_to_ref_aligned.bam.bai

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai
        
        tabix -p vcf clinvar.vcf.gz

        /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref genome_reference.fasta --reads ~{file_label}_consensus_to_ref_aligned.bam --output_vcf  ~{file_label}_consensus_to_ref_dv_raw_variants.vcf.gz
        
        /opt/deepvariant/bin/vcf_stats_report --input_vcf ~{file_label}_consensus_to_ref_dv_raw_variants.vcf.gz --outfile_base ~{file_label}
        
        bcftools annotate -c ID,INFO -a clinvar.vcf.gz ~{file_label}_consensus_to_ref_dv_raw_variants.vcf.gz | bcftools csq -f genome_reference.fasta -g ~{gff} | bcftools reheader -s temp_header.txt | bcftools view -Ov -o ~{file_label}_consensus_to_ref_raw_variants_annotated.vcf.gz

        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ]\n" ~{file_label}_consensus_to_ref_raw_variants_annotated.vcf.gz > ~{file_label}_dv_raw_variant_summary.tsv

    >>>

    output {
        Array[File] raw_vcf = glob("*.vcf.gz")
        Array[File] raw_vcf_report = glob("*.html")
    }

    runtime {
        docker: "google/deepvariant:~{deepvariant_version}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
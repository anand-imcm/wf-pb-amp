version 1.0

# align clustered hifi reads to reference and generate the final bam using pbaa bampaint
task HifiReadsVarCall {
    
    input {
        File raw_hifi_to_reference_alignment_bam
        File raw_hifi_to_reference_alignment_index
        File genome_reference
        File clinvar
        File gff
        String file_label
        String deepvariant_version = "1.5.0"
        Int deepvariant_num_shards = 8
    }  

    command <<<
        set -euo pipefail

        ln -s ~{genome_reference} genome_reference.fasta

        ln -s ~{clinvar} clinvar.vcf.gz

        ln -s ~{raw_hifi_to_reference_alignment_bam} ~{file_label}_raw_hifi_to_reference_alignment.bam

        ln -s ~{raw_hifi_to_reference_alignment_index} ~{file_label}_raw_hifi_to_reference_alignment.bam.bai

        samtools faidx genome_reference.fasta -o genome_reference.fasta.fai
        
        tabix -p vcf clinvar.vcf.gz

        /opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --num_shards ~{deepvariant_num_shards} \
            --ref genome_reference.fasta \
            --reads ~{file_label}_raw_hifi_to_reference_alignment.bam \
            --output_vcf  ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz
        
        /opt/deepvariant/bin/vcf_stats_report \
            --input_vcf ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz \
            --outfile_base ~{file_label}
        
        bcftools annotate -c ID,INFO -a clinvar.vcf.gz ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz | bcftools csq -f genome_reference.fasta -g ~{gff} | bcftools view -Ov -o ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz

        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ]\n" ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz > ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated_summary.tsv

    >>>

    output {
        File raw_hifi_to_reference_alignment_all_variants_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_annotated_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_stats = file_label + ".visual_report.html"
        File raw_hifi_to_reference_alignment_all_variants_annotated_summary = file_label + "_raw_hifi_to_reference_alignment_all_variants_annotated_summary.tsv"
    }

    runtime {
        docker: "google/deepvariant:~{deepvariant_version}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
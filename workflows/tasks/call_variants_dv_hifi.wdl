version 1.0

# call variants using Google deep variant
task HifiReadsVarCallDV {
    
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
        
        bcftools annotate -c ID,INFO -a clinvar.vcf.gz ~{file_label}_raw_hifi_to_reference_alignment_all_variants.vcf.gz | bcftools csq -p a -f genome_reference.fasta -g ~{gff} | bcftools view -Oz -o ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz

        tabix -p vcf ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz

        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ\t%GT:%GQ:%DP:%AD:%VAF:%PL:%BCSQ]\n" ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz > ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated_summary.tsv

        modified_header_all_vars=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
        
        # replacing the header infile
        sed -i "1s/.*/$modified_header_all_vars/" ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated_summary.tsv

        bcftools view -f PASS ~{file_label}_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz -Oz -o ~{file_label}_raw_hifi_to_reference_alignment_pass_variants_annotated.vcf.gz

        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ\t%GT:%GQ:%DP:%AD:%VAF:%PL:%BCSQ]\n" ~{file_label}_raw_hifi_to_reference_alignment_pass_variants_annotated.vcf.gz > ~{file_label}_raw_hifi_to_reference_alignment_pass_variants_annotated_summary.tsv
        
        modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_pass_variants_annotated_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
        
        # replacing the header infile
        sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_pass_variants_annotated_summary.tsv
    >>>

    output {
        File raw_hifi_to_reference_alignment_all_variants_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_annotated_vcf = file_label + "_raw_hifi_to_reference_alignment_all_variants_annotated.vcf.gz"
        File raw_hifi_to_reference_alignment_all_variants_annotated_summary = file_label + "_raw_hifi_to_reference_alignment_all_variants_annotated_summary.tsv"
        File raw_hifi_to_reference_alignment_pass_variants_annotated_vcf = file_label + "_raw_hifi_to_reference_alignment_pass_variants_annotated.vcf.gz"
        File raw_hifi_to_reference_alignment_pass_variants_annotated_summary = file_label + "_raw_hifi_to_reference_alignment_pass_variants_annotated_summary.tsv"
        File raw_hifi_to_reference_alignment_all_variants_stats = file_label + ".visual_report.html"
    }

    runtime {
        docker: "google/deepvariant:~{deepvariant_version}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
version 1.0

# get ontarget varaints from the Google deep variant call set
task HifiOnTargetVarsDV {
    
    input {
        File raw_hifi_to_reference_alignment_pass_variants_annotated_vcf
        File bed
        String file_label
        String docker        
    }  

    command <<<
        set -euo pipefail

        ln -s ~{raw_hifi_to_reference_alignment_pass_variants_annotated_vcf} pass.annotated.vcf.gz

        tabix -p vcf pass.annotated.vcf.gz

        bedtools intersect -header -a pass.annotated.vcf.gz -b ~{bed} -wa | bgzip -c > ~{file_label}_raw_hifi_to_reference_alignment_ontarget_pass_variants.vcf.gz
        
        tabix -p vcf ~{file_label}_raw_hifi_to_reference_alignment_ontarget_pass_variants.vcf.gz
        
        bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ\t%GT:%GQ:%DP:%AD:%VAF:%PL:%BCSQ]\n" ~{file_label}_raw_hifi_to_reference_alignment_ontarget_pass_variants.vcf.gz > ~{file_label}_raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary.tsv

        modified_header=$(head -n1 ~{file_label}_raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary.tsv | sed 's/\[[0-9]*\]//g; s/#//')
        
        # replacing the header infile
        sed -i "1s/.*/$modified_header/" ~{file_label}_raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary.tsv

    >>>

    output {
        File raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_vcf = file_label + "_raw_hifi_to_reference_alignment_ontarget_pass_variants.vcf.gz"
        File raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary = file_label +"_raw_hifi_to_reference_alignment_ontarget_pass_variants_annotated_summary.tsv"

    }

    runtime {
        docker: "~{docker}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
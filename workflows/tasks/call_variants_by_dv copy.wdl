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
        
        # python /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/scripts/variant_call_deepvariant.py \
        #     --inbam ~{file_label}_consensus_to_ref_aligned.bam \
        #     --ref genome_reference.fasta \
        #     --clinvar clinvar.vcf.gz \
        #     --gff ~{gff} \
        #     --target ~{bed} \
        #     --prefix ~{file_label}

        python <<CODE
        #! python3

        import os, pysam, argparse
        import pandas as pd

        tsvs_all_variants=[]
        tsvs_all_variants_summary = f'~{file_label}_variant_summary.tsv'

        with pysam.AlignmentFile(~{file_label}_consensus_to_ref_aligned.bam) as inbam:
            if (inbam.check_index):
                if (inbam.mapped > 0) :
                    for rec in inbam :
                        obam = f'{rec.query_name}_reference.bam'
                        vcf = f'{rec.query_name}_reference_raw.vcf.gz'
                        vcf_annotated = f'{rec.query_name}_reference_annotated.vcf.gz'
                        variants_tsv = f'{rec.query_name}_reference_variants.tsv'
                        variants_ontarget = f'{rec.query_name}_reference_ontarget_variants.vcf.gz'
                        variants_ontarget_tsv = f'{rec.query_name}_reference_ontarget_variants.tsv'
                        pysam.AlignmentFile(obam, 'wb', header=inbam.header).write(rec)
                        pysam.index(obam)
                        with open('temp_header.txt', 'w') as temp_header_file:
                            temp_header_file.write(rec.query_name)
                        os.system(f'/opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref genome_reference.fasta --reads {obam} --output_vcf {vcf}')
                        os.system(f'/opt/deepvariant/bin/vcf_stats_report --input_vcf {vcf} --outfile_base {rec.query_name}')
                        os.system(f'tabix -p vcf {vcf}')
                        os.system(f'bcftools annotate -c ID,INFO -a clinvar.vcf.gz {vcf} | bcftools csq -f genome_reference.fasta -g ~{gff} | bcftools reheader -s temp_header.txt | bcftools view -Ov -o {vcf_annotated}')
                        os.system(f'bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ]\n" {vcf_annotated} > {variants_tsv}')
                        tsvs_all_variants.append(variants_tsv)
                else :
                    print(f"Skipping: ~{file_label}. The total number of mapped alignments is 0")

        patt = 'sample-(?P<Sample>.*)_guide-(?P<Target>.*)_cluster-(?P<Cluster>[0-9]+)_ReadCount-(?P<Numreads>[0-9]+)'

        def read_tsv( tsv ):
            tbl = pd.read_csv( tsv, sep='\t' )
            if tbl:
                tbl.columns = tbl.columns.str.replace( '^.*\\]', '', regex=True ).str.split(':').str[-1]
                if tbl.empty:
                    tbl.loc[0] = {'SAMPLE' : tsv, 'CHROM': 'NoVariants'}
                
                #fill in columns with missing cluster ("sample") name
                tbl.SAMPLE = tbl.SAMPLE[tbl.SAMPLE.notnull()].iloc[0]
                clusterInfo = tbl.SAMPLE.str.extract(patt)
                res = pd.concat([tbl, clusterInfo], axis=1).drop(columns='SAMPLE').fillna('.')
                res = res.reindex([ 'Sample','Target','Cluster','Numreads','CHROM','POS','ID','REF','ALT','BCSQ'], axis=1)
                return res

        if tsvs_all_variants:
            pd.concat(map(read_tsv,tsvs_all_variants), axis=0)\
            .fillna('.')\
            .to_csv(tsvs_all_variants_summary, sep="\t", index=False)

        else:
            pd.DataFrame(columns=['Sample','Target','Cluster','Numreads','CHROM','POS','ID','REF','ALT','BCSQ'])\
            .fillna('.')\
            .to_csv(tsvs_all_variants_summary, sep="\t", index=False)

        CODE
        
    >>>

    output {
        File variant_summary = file_label + "_variant_summary.tsv"
        Array[File] annotated_vcf = glob("*_annotated.vcf.gz")
        Array[File] raw_vcf = glob("*_raw.vcf.gz")
        Array[File] raw_vcf_report = glob("*.html")
    }

    runtime {
        docker: "google/deepvariant:~{deepvariant_version}"
        memory: "32G"
        disks: "local-disk 30 HDD"
    }
}
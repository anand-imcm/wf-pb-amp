#! python3

import os, pysam, argparse
import pandas as pd

parser = argparse.ArgumentParser(prog='variant_call.py', description='Call variants per concensus')
parser.add_argument('-i','--inbam',type=str, help='pbaa consensus to reference aligned bam', required=True)
parser.add_argument('-r','--ref',type=str, help='Genome reference fasta (must be indexed)', required=True)
parser.add_argument('-c','--clinvar',type=str, help='Clinvar vcf', required=True)
parser.add_argument('-g','--gff',type=str, help='Consequence gff', required=True)
parser.add_argument('-t','--target',type=str, help='Target region coordinates in bed format', required=True)
parser.add_argument('-p','--prefix',type=str, help='Output variant summary from all clusters in tsv format', required=True)

args = parser.parse_args()

tsvs_all_variants=[]
tsvs_on_target=[]
tsvs_all_variants_summary = f'{args.prefix}_variant_summary.tsv'
tsvs_on_target_summary = f'{args.prefix}_variant_on_target_summary.tsv'

with pysam.AlignmentFile(args.inbam) as inbam:
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
                os.system(f'bcftools mpileup -f {args.ref} {obam} | bcftools call -mv -Oz -o {vcf}')
                os.system(f'tabix -p vcf {vcf}')
                os.system(f'bcftools annotate -c ID,INFO -a {args.clinvar} {vcf} | bcftools csq -f {args.ref} -g {args.gff} | bcftools reheader -s temp_header.txt | bcftools view -Ov -o {vcf_annotated}')
                os.system(f'bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ]\n" {vcf_annotated} > {variants_tsv}')
                os.system(f'bedtools intersect -header -a {vcf_annotated} -b {args.target} -wa | bgzip -c > {variants_ontarget}')
                os.system(f'tabix -p vcf {variants_ontarget}')
                os.system(f'bcftools query -Hu -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SAMPLE\t%INFO/BCSQ]\n" {variants_ontarget} > {variants_ontarget_tsv}')
                tsvs_all_variants.append(variants_tsv)
                tsvs_on_target.append(variants_ontarget_tsv)
        else :
            print(f"Skipping: {args.prefix}. The total number of mapped alignments is 0")


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
    pd.concat( map( read_tsv, tsvs_all_variants ), axis=0 )\
    .fillna('.')\
    .to_csv(tsvs_all_variants_summary, sep="\t", index=False)

    pd.concat( map( read_tsv, tsvs_on_target ), axis=0 )\
    .fillna('.')\
    .to_csv(tsvs_on_target_summary, sep="\t", index=False)
else:
    pd.DataFrame(columns=['Sample','Target','Cluster','Numreads','CHROM','POS','ID','REF','ALT','BCSQ'])\
    .fillna('.')\
    .to_csv(tsvs_all_variants_summary, sep="\t", index=False)

    pd.DataFrame(columns=['Sample','Target','Cluster','Numreads','CHROM','POS','ID','REF','ALT','BCSQ'])\
    .fillna('.')\
    .to_csv(tsvs_on_target_summary, sep="\t", index=False)
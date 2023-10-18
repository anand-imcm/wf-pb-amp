#! python3

import pysam, argparse

parser = argparse.ArgumentParser(prog='create_guide_fasta.py', description='Create guide reference fasta')
parser.add_argument('-r','--ref',type=str, help='Genome reference fasta', required=True)
parser.add_argument('-t','--target',type=str, help='Target region coordinates in bed format', required=True)
parser.add_argument('-o','--out',type=str, help='Output guide reference fasta file name', required=True)

args = parser.parse_args()

# usage: python scripts/create_guide_fasta.py --ref /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/data/new-set/reference/Homo_sapiens.GRCh38.release110.dna.chromosome.1.fa --target /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/data/new-set/reference/targets.bed --out /home/anand/Documents/aspire-files/data-oxford/terra.bio/wf-pb-amp/data/new-set/reference/gba1.gbap1.guide.from.Homo_sapiens.GRCh38.dna.chromosome.1.fasta

regions={}
with open(args.target) as bed:
   for coord in bed:
        rec = coord.strip().split("\t")
        regions[rec[3]]=f'{rec[0]}:{rec[1]}-{rec[2]}'

with pysam.FastaFile(args.ref) as infa, open(args.out, 'w') as outfa:
    for label,region in regions.items():
        sequence = infa.fetch( region=region )
        outfa.write( f'>{label}\n{sequence}\n' )
pysam.faidx( args.out )
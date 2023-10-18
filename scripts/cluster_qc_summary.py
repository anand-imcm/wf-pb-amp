#! python3

import pandas as pd
import argparse

parser = argparse.ArgumentParser(prog='cluster_qc_summary.py', description='Create summary using the clustered reads')
parser.add_argument('-r','--readinfo',type=str, help='Genome reference fasta', required=True)
parser.add_argument('-d','--demuxreport',type=str, help='Target region coordinates in bed format', required=True)
parser.add_argument('-p','--prefix',type=str, help='sample name or summary file prefix', required=True)

args = parser.parse_args()

pbaa_info = args.readinfo
demux_report = args.demuxreport
output = f'{args.prefix}_clusterQC_report.tsv'

# Load cluster_info from pbaa_info
cluster_info = pd.read_csv(pbaa_info, sep='\\s', engine='python', usecols=[0, 1, 9], names=['read', 'target', 'cluster'], index_col=0)
cluster_info.index = cluster_info.index.str.rsplit('/', n=1).str[0]

# Load end_calls data from demux_report
end_calls = pd.read_csv(demux_report, sep='\t', usecols=[0, 7, 35, 36], names=['read', 'qual', 'fwd', 'rev'], index_col=0)

# Join with cluster_info
end_calls = end_calls.join(cluster_info)

# Group and calculate statistics
counts = end_calls.groupby(['target', 'cluster']).describe(percentiles=[0.9]).sort_index()

# Calculate frequencies
freqs = counts[('qual', 'count')] / len(end_calls)
freqs.name = ('', 'frequency')
counts = counts.join(freqs).set_index(('','frequency'), append=True)
counts.index.names = counts.index.names[:-1] + ['frequency']

# Filter and save the result
counts.query('frequency >= 0.02').to_csv(output, sep='\t', float_format='%.3f')
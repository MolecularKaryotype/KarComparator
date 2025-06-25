import Analyses_UTILS as AU
import pandas as pd
import os
import argparse
import sys


## for log output
class Tee(object):
    def __init__(self, *streams):
        self.streams = streams

    def write(self, data):
        for s in self.streams:
            s.write(data)
            s.flush()  # Ensure it's written immediately

    def flush(self):
        for s in self.streams:
            s.flush()

parser = argparse.ArgumentParser()
parser.add_argument('--matching_allowance', type=int, default=200000, help='set matching allowance for KarCheck')
parser.add_argument('--approx_allowance', type=int, default=50000, help='set indel minimum size for KarCheck')
parser.add_argument('--karsim_kt_dir', type=str, help='path to DIR of simulation karyotypes')
parser.add_argument('--karsim_edgelabel_dir', type=str, help="path to DIR of simulation discordant edge labels")
parser.add_argument('--metadata_dir', type=str, help="path to DIR containing metadata such as reference labeling")
parser.add_argument('--output_dir', type=str, help="path to this program's output")
parser.add_argument('--clusters_dir', type=str, help="path to DIR containing the independent clusters")
parser.add_argument('--omkar_output_dir', type=str, help="path to DIR containing omkar outputs")
parser.add_argument('--simulation_bionano_output', type=str, help="path to DIR containing bionano output from the simulated TRUTH karyotypes")
args = parser.parse_args()


## running KarCheck
karcheck_matching_allowance = args.matching_allowance
karcheck_approx_allowance = args.approx_allowance
karsim_data_dir = args.karsim_kt_dir
karsim_history_edge_dir = args.karsim_edgelabel_dir
metadata_dir = args.metadata_dir
output_dir = args.output_dir
os.makedirs(os.path.join(output_dir, 'analyses_summary/'), exist_ok=True)

AU.data_folder = args.clusters_dir
AU.omkar_log_dir = args.omkar_output_dir
AU.karsim_file_prefix = karsim_data_dir
AU.karsim_history_edges_folder = karsim_history_edge_dir
AU.forbidden_region_file = f"{metadata_dir}/acrocentric_telo_cen.bed"
AU.bionano_data_folder = args.simulation_bionano_output
cn_bin_file = f"{metadata_dir}/cn_bins_50kbp.txt"
karcheck_log_file = f"{output_dir}/analyses_summary/summary_stats.txt"

## writing streams to both stdout and log file
log_write = open(karcheck_log_file, 'w')
tee = Tee(sys.stdout, log_write)
sys.stdout = tee

df = AU.prep_df()
df, s1 = AU.KarCheck(df, approximate_allowance=karcheck_approx_allowance, matching_allowance=karcheck_matching_allowance)
print(s1)
df_nonevent = df[~df['is_event_cluster']]
df_event = df[df['is_event_cluster']]

print('\nSummary stats on number of clusters')
print('Cases with wrong number of chromosome reconstructed: ')
print(f"number of event clusters: {len(df_event)}")
print(f"number of non-event clusters: {len(df_nonevent)}")
print(f"number of non-event cluster with SV reconstructed (FP-only cluster): {(df_nonevent['nT_FP_edges'] != 0).sum()}")
AU.summary_statistics_by_complexity_group(df)

total_nt = 0
for row_idx, df_row in df.iterrows():
    for event_id, status in df_row['event_status'].items():
        total_nt += 1
print(f'\ntotal nT events simulated: {total_nt}')

print('\n--------------------------------------------------------------------------------')
print('## Recall by SV-edge simulated, separated by SV types')
AU.summary_statistics_by_event_types(df)

### reading bed file + KarCheck output table
### for tallying the cnv-missed partially-supported event calls
### !not incorporated! into OMKar output
df['bed_file_residual_edge'] = df.apply(lambda row: AU.rowwise_match_residual_edge_in_bed_file(row), axis=1)
print(f"Number of cnv-missed partially-supported SVs in BED file, not in OMKar reconstruction: {df['bed_file_residual_edge'].apply(len).sum()}")

print('\nnT report, event type recall <= 6', end="")
AU.summary_statistics_by_event_types(df[df['nT_complexity'] <= 6])
print('\nnT report, event type recall > 6', end="")
AU.summary_statistics_by_event_types(df[df['nT_complexity'] > 6])
print('\nT report, event type recall <= 6', end="")
AU.summary_statistics_by_event_types(df[df['TnT_complexity'] <= 6], status_column='T_event_status')
print('\nT report, event type recall > 6', end="")
AU.summary_statistics_by_event_types(df[df['TnT_complexity'] > 6], status_column='T_event_status')
print('\nT report, event type recall ALL', end="")
AU.summary_statistics_by_event_types(df, status_column='T_event_status')

print('\n--------------------------------------------------------------------------------')
print('## Number of Chromosome reconstructed analyses')

print('\nCases with wrong number of chromosome reconstructed: ')
AU.report_aneuploidy_anomaly_cases(df, 'nT_aneuploidy_comparison')

print('\nnT-only event aneuploidy')
AU.total_cluster_with_wrong_number_of_chrom(df_event[df_event['T_event_status'].apply(len) == 0], 'nT_aneuploidy_comparison')
print('\nall nonevent aneuploidy')
AU.total_cluster_with_wrong_number_of_chrom(df_nonevent, 'TnT_aneuploidy_comparison')

print('\nall event aneuploidy')
terminal_df_with_nonterminal_event = df_event[(df_event['T_event_status'].apply(lambda d: len(d) != 0)) & (df_event['nT_event_status'].apply(lambda d: len(d) != 0))]
AU.total_cluster_with_wrong_number_of_chrom(terminal_df_with_nonterminal_event, 'TnT_aneuploidy_comparison')

print('\nall nT aneuploidy')
df_with_nonterminal_event = df_event[df_event['nT_event_status'].apply(lambda d: len(d) != 0)]
AU.total_cluster_with_wrong_number_of_chrom(df_with_nonterminal_event, 'TnT_aneuploidy_comparison')

print('\nT-only aneuploidy')
df_with_nonterminal_event_without_terminal = df_event[(df_event['nT_event_status'].apply(lambda d: len(d) != 0)) & (df_event['T_event_status'].apply(lambda d: len(d) == 0))]
AU.total_cluster_with_wrong_number_of_chrom(df_with_nonterminal_event_without_terminal, 'TnT_aneuploidy_comparison')

print('\nall T aneuploidy')
AU.total_cluster_with_wrong_number_of_chrom(df[df['T_event_status'].apply(len) != 0], 'TnT_aneuploidy_comparison')

df.to_csv(f"{output_dir}/analyses_summary/whole_analysis_table.csv", index=False)

## Stable CN comparison, disabled for better runtimes
print('\n--------------------------------------------------------------------------------')
print('## CN Concordance Analyses')
print('CN started')
df[['karsim_CN', 'omkar_CN']] \
    = df.apply(lambda row: pd.Series(AU.iterative_get_cn_with_bins(row, cn_file_name=cn_bin_file)), axis=1)
df['CN_jaccard'] = df.apply(lambda row: AU.rowwise_cn_jaccard_sim_by_cluster(row, rounding_allowance=0.05), axis=1)
print(f"CN Jaccard Score: {df['CN_jaccard'].mean()}")
df.to_csv(f"{output_dir}/analyses_summary/whole_analysis_table_with_cn_005.csv", index=False)

## restore stdout writing stream
sys.stdout = sys.__stdout__
log_write.close()


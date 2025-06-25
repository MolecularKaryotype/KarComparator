import argparse
import pandas as pd
import os

os.makedirs('comparisons', exist_ok=True)

parser = argparse.ArgumentParser()
parser.add_argument('run1_version', type=str, help='version name of run1')
parser.add_argument('run2_version', type=str, help='version name of run2')
args = parser.parse_args()

run1_version = args.run1_version
run2_version = args.run2_version

run1_csv_path = 'builds/{}/analyses_summary/whole_analyses_table.csv'.format(run1_version)
run2_csv_path = 'builds/{}/analyses_summary/whole_analyses_table.csv'.format(run2_version)

df1 = pd.read_csv(run1_csv_path)
df2 = pd.read_csv(run2_csv_path)
df1 = df1[['file_name', 'cluster', 'origin_chr', 'n_path_diff', 'SV_missed']]
df2 = df2[['file_name', 'cluster', 'origin_chr', 'n_path_diff', 'SV_missed']]

merged_df = pd.merge(df1, df2, on=['file_name', 'cluster'], suffixes=('_df1', '_df2'))
merged_df = merged_df[['file_name', 'cluster', 'origin_chr_df1', 'origin_chr_df2', 'n_path_diff_df1', 'n_path_diff_df2', 'SV_missed_df1', 'SV_missed_df2']]
merged_df['SV_missed_delta'] = merged_df['SV_missed_df2'] - merged_df['SV_missed_df1']

path_diff_df = merged_df[merged_df['n_path_diff_df1'] != merged_df['n_path_diff_df2']]
path_diff_df = path_diff_df.sort_values(by=['n_path_diff_df1'])
path_diff_df = path_diff_df.sort_values(by=['n_path_diff_df2'])
path_diff_df = path_diff_df.reindex(path_diff_df['n_path_diff_df2'].abs().sort_values(ascending=False).index)
path_diff_df.to_csv('comparisons/{}vs{}_path_diffs.csv'.format(run1_version, run2_version), index=False)

sv_missed_df = merged_df[merged_df['SV_missed_df1'] != merged_df['SV_missed_df2']]
sv_missed_df = sv_missed_df.sort_values(by=['SV_missed_delta'], ascending=False)
sv_missed_df.to_csv('comparisons/{}vs{}_SV_missed.csv'.format(run1_version, run2_version), index=False)

from COMPARISON_with_graphs import *
from read_KarSimulator_output import *
from KarUtils import *
from bionano_output_parser import *

import os
import pandas as pd
import numpy as np
import re
"""
README: USAGE PIPELINE
[Fill up the variables below]
prep_df(): file_name, cluster, n_origin_chr, origin_chr, n_path_karsim, n_path_omkar, event_chr, histories, n_events, is_event_cluster
KarCheck(df): 
"""


### MUST be overwritten during IMPORT
data_folder = ''
omkar_log_dir = ''
karsim_file_prefix = ''
karsim_history_edges_folder = ''
forbidden_region_file = ''
bionano_data_folder = ''

SIG_DIST_THRESHOLD = 100000


def form_graph_from_cluster_file(input_df_row):
    file_prefix = data_folder
    file_suffix = '.txt'
    cluster_file_name = file_prefix + input_df_row['file_name'] + 'cluster_' + str(input_df_row['cluster']) + file_suffix
    graph = form_graph_from_cluster(cluster_file_name)
    return graph


def iterative_count_events(df_row):
    event_dict = df_row['histories']
    tot = 0
    for key, value in event_dict.items():
        tot += value
    return tot


def prep_df():
    """
    populate the following columns:
    file_name, cluster, n_origin_chr, origin_chr, n_path_karsim, n_path_omkar, event_chr, histories, n_events, is_event_cluster
    :return:
    """
    files = []
    for file in os.listdir(data_folder):
        files.append(file)
    files.sort()

    ## Extract basic cluster file's info
    data = []
    for file in files:
        file_name = file.split('cluster')[0]
        with open(data_folder + file) as fp_read:
            line1 = fp_read.readline()
            matches = re.findall(r'<(.*?)>', line1)
            origins = eval(matches[1])
            new_data = {'file_name': file_name,
                        'cluster': matches[0],
                        'n_origin_chr': len(origins),
                        'origin_chr': origins,
                        'n_path_karsim': int(matches[2]),
                        'n_path_omkar': int(matches[3])}
        data.append(new_data)

    df = pd.DataFrame(data)

    ## n_Chr difference
    df['n_path_diff'] = df['n_path_omkar'] - df['n_path_karsim']

    ## get event_chr
    df['event_chr'] = df['file_name'].apply(lambda x: list(get_event_chr(karsim_file_prefix + x + '.kt.txt')))
    df['event_chr'] = df['event_chr'].apply(lambda x: [entry[:-1] for entry in x])

    df['histories'] = df.apply(lambda x: get_history_events(karsim_file_prefix + x['file_name'] + '.kt.txt',
                                                            x['origin_chr']), axis=1)
    df['n_events'] = df.apply(lambda row: iterative_count_events(row), axis=1)
    df['is_event_cluster'] = df.apply(lambda row: any(origin_chr in row['event_chr'] for origin_chr in row['origin_chr']), axis=1)

    return df


def KarCheck(df, approximate_allowance=50000, matching_allowance=200000, karsim_terminal_distance=5000, omkar_terminal_distance=200000):
    """
    This is aimed to be the finalized version for KarCheck
    require running prep_df() before running this
    :param df:
    :return:
    """
    def rowwise_label_edges(df_row):
        event_sv_edges = case_event_edges[df_row['file_name']]
        df_row['graph'].add_sv_label_to_karsim_edges(event_sv_edges)

    def rowwise_get_terminal_event_id(df_row):
        graph = df_row['graph']
        karsim_T, karsim_nT, omkar_T, omkar_nT = graph.get_terminally_labeled_residual_transitions(karsim_terminal_distance, omkar_terminal_distance)
        terminal_event_ids = []
        for entry in karsim_T:
            for edge, ids in entry.items():
                terminal_event_ids += ids
        return list(set(terminal_event_ids))

    def rowwise_get_nonterminal_omkar_residuals(df_row):
        graph = df_row['graph']
        karsim_T, karsim_nT, omkar_T, omkar_nT = graph.get_terminally_labeled_residual_transitions(karsim_terminal_distance, omkar_terminal_distance)
        return len(omkar_nT)

    def rowwise_read_smap(df_row):
        smap_file_path = f"{bionano_data_folder}/{df_row['file_name']}/exp_refineFinal1_merged_filter_inversions.smap"
        return smap_to_df(smap_file_path)

    def rowwise_read_bed(df_row):
        bed_file_path = f"{omkar_log_dir}/{df_row['file_name']}/{df_row['file_name']}_SV.bed"
        return read_bed_file(bed_file_path)

    def rowwise_check_edge_in_smap(df_row, distance, status_col):
        def rename_chr(input_chr):
            new_chr = input_chr.upper().replace('CHR', '')
            if new_chr == 'X':
                return 23
            elif new_chr == 'Y':
                return 24
            else:
                return int(new_chr)

        caught_event_status = df_row[status_col]
        smap_df = df_row['SMAP']
        smap_status = {}
        # there can be multiple SV-events per cluster
        for event_id, status in caught_event_status.items():
            karsim_edge_in_smap = []
            omkar_edge_in_smap = []
            for edge in status['karsim_edge']:
                # there can be multiple edges per SV-event
                filtered_df = filter_chroms(smap_df, rename_chr(edge[0]), rename_chr(edge[2]))
                if 'inversion' not in event_id:
                    filtered_df = approx_edge_search(filtered_df, edge[1], edge[3], distance=distance)
                    if not filtered_df.empty:
                        karsim_edge_in_smap.append(True)
                    else:
                        karsim_edge_in_smap.append(False)
                else:
                    # for inversion edges, the two positions may be separated into two partial edges (or two full edges each with a confident position)
                    filtered_df_pos1, filtered_df_pos2 = approx_edge_search_inv(filtered_df, edge[1], edge[3], distance)
                    if (not filtered_df_pos1.empty) and (not filtered_df_pos2.empty):
                        karsim_edge_in_smap.append(True)
                    else:
                        karsim_edge_in_smap.append(False)
            for edge in status['omkar_edge']:
                filtered_df = filter_chroms(smap_df, rename_chr(edge[0]), rename_chr(edge[2]))
                if 'inversion' not in event_id:
                    filtered_df = approx_edge_search(filtered_df, edge[1], edge[3], distance=distance)
                    if not filtered_df.empty:
                        omkar_edge_in_smap.append(True)
                    else:
                        omkar_edge_in_smap.append(False)
                else:
                    # for inversion edges, the two positions may be separated into two partial edges (or two full edges each with a confident position)
                    filtered_df_pos1, filtered_df_pos2 = approx_edge_search_inv(filtered_df, edge[1], edge[3], distance)
                    if (not filtered_df_pos1.empty) and (not filtered_df_pos2.empty):
                        omkar_edge_in_smap.append(True)
                    else:
                        omkar_edge_in_smap.append(False)
            smap_status[event_id] = {'karsim_edge': karsim_edge_in_smap, 'omkar_edge': omkar_edge_in_smap}
        return smap_status

    def summarize_aneuploidy_stats(df_row, aggregated_stats):
        aneuploidy_info = df_row['nT_aneuploidy_comparison']
        for categority, count in aneuploidy_info.items():
            aggregated_stats[categority] += count

    case_event_edges = {}
    for case in df['file_name'].unique():
        sv_label_filepath = f"{karsim_history_edges_folder}/{case}.history_sv.txt"
        case_event_edges[case] = read_history_edges_intermediate_file(sv_label_filepath)
    df['graph'] = df.apply(form_graph_from_cluster_file, axis=1)
    df['unmodified_graph'] = df.apply(form_graph_from_cluster_file, axis=1)
    df['graph'].apply(lambda g: g.remove_all_segment_edges())
    df['graph'].apply(lambda g: g.remove_approximate_transition_edges(approximate_allowance))
    ## label all karsim edges and categorize events into terminal vs nonterminal (T vs nT)
    df.apply(rowwise_label_edges, axis=1)
    df['T_event_ids'] = df.apply(lambda row: rowwise_get_terminal_event_id(row), axis=1)
    ## create event status object to track all events' edges matching status and distances
    df['event_status'] = df['graph'].apply(lambda g: tally_graph_edge_labels(g))
    df['event_status'].apply(lambda l: label_status_log_initial_count(l))
    df['event_status'].apply(lambda l: label_status_populate_edge_fields(l))
    ## matching
    df['same_edge_status'] = df['graph'].apply(lambda g: g.prune_same_edges())
    df.apply(lambda row: update_label_status(row['event_status'], row['same_edge_status']), axis=1)
    df['matched_edge_status'] = df['graph'].apply(lambda g: g.match_transition_edges(matching_allowance))
    df.apply(lambda row: update_label_status(row['event_status'], row['matched_edge_status']), axis=1)
    df[['karsim_residual', 'omkar_residual']] = df['graph'].apply(lambda g: pd.Series(g.get_missed_transition_edges()))
    ## split event status into T vs nT
    df[['T_event_status', 'nT_event_status']] = df.apply(lambda row: pd.Series(split_event_status_into_terminal_nonterminal(row['T_event_ids'], row['event_status'])), axis=1)
    ## extract events that are completely caught
    df[['nT_fully_caught_event_status', 'nT_not_fully_caught_event_status']] = df['nT_event_status'].apply(lambda l: pd.Series(split_event_status_into_caught_uncaught_events(l)))
    ## summary statistics based on edges and events (nT only)
    df['nT_complexity'] = df['nT_event_status'].apply(lambda l: complexity_by_event_status(l))
    df['TnT_complexity'] = df['event_status'].apply(lambda l: complexity_by_event_status(l))
    compare_aneuploidy_for_nT_clusters(df)
    compare_aneuploidy_for_TnT_clusters(df)
    df[['nT_TP_events', 'nT_FN_events', 'nT_TP_edges', 'nT_FN_edges']] = df['nT_event_status'].apply(lambda l: pd.Series(event_status_get_caught_uncaught_counts(l)))
    df['nT_FP_edges'] = df.apply(lambda row: rowwise_get_nonterminal_omkar_residuals(row), axis=1)
    ## compute final scores for combined SV edges and aneuploidy results
    compute_row_recall_precision_jaccard(df, 'nT',
                                         TP=df['nT_TP_edges'] + df['nT_aneuploidy_comparison'].apply(lambda x: x['TP']),
                                         FP=df['nT_FP_edges'] + df['nT_aneuploidy_comparison'].apply(lambda x: x['FP']),
                                         FN=df['nT_FN_edges'] + df['nT_aneuploidy_comparison'].apply(lambda x: x['FN']))
    output_str = "SV-edge summary stats:\n"
    output_str += f"TP (total): {(df['nT_TP_edges'] + df['nT_aneuploidy_comparison'].apply(lambda x: x['TP'])).sum()}\n"
    output_str += f"FP (total): {(df['nT_FP_edges'] + df['nT_aneuploidy_comparison'].apply(lambda x: x['FP'])).sum()}\n"
    output_str += f"FN (total): {(df['nT_FN_edges'] + df['nT_aneuploidy_comparison'].apply(lambda x: x['FN'])).sum()}\n"
    output_str += f"Avg Jaccard Case: {df[df['nT_jaccard'] != -1]['nT_jaccard'].mean()}\n"
    output_str += f"Avg Recall Case: {df[df['nT_recall'] != -1]['nT_recall'].mean()}\n"
    output_str += f"Avg Precision Case: {df[df['nT_precision'] != -1]['nT_precision'].mean()}\n"
    output_str += f"total number of clusters: {df.shape[0]}\n"

    ## aneuploidy summaries (already integrated into scoring for each cluster)
    aneuploiody_stats = {'TN': 0, 'TP': 0, 'FN': 0, 'FP': 0}
    df.apply(lambda row: summarize_aneuploidy_stats(row, aneuploiody_stats), axis=1)
    output_str += f"aneuploidy stats: {aneuploiody_stats}"

    ## analyses of matching distances
    df['SMAP'] = df.apply(lambda row: rowwise_read_smap(row), axis=1)
    df['BED'] = df.apply(lambda row: rowwise_read_bed(row), axis=1)
    df['edges_in_smap'] = df.apply(lambda row: rowwise_check_edge_in_smap(row, distance=50000, status_col='nT_fully_caught_event_status'), axis=1)
    df['T_edges_in_smap'] = df.apply(lambda row: rowwise_check_edge_in_smap(row, distance=50000, status_col='T_event_status'), axis=1)

    return df, output_str


def report_aneuploidy_anomaly_cases(df, col):
    def rowwise_report_aneuploidy_anomaly(df_row):
        FP = df_row[col]['FP']
        FN = df_row[col]['FN']
        if FP + FN != 0:
            print(f"case: {df_row['file_name']}; cluster: {df_row['cluster']}; origin_chr: {df_row['origin_chr']}; FP: {FP}; FN: {FN}")
    df.apply(lambda row: rowwise_report_aneuploidy_anomaly(row), axis=1)


def summary_statistics_by_complexity_group(df):
    jaccard_df = df[df['nT_jaccard'] != -1]
    precision_df = df[df['nT_precision'] != -1]
    recall_df = df[df['nT_recall'] != -1]

    def output_score_for_each_complexity(sub_df, score_col):
        low_df = sub_df[(sub_df['nT_complexity'] >= 1) & (sub_df['nT_complexity'] <= 2)]
        mid_df = sub_df[(sub_df['nT_complexity'] >= 3) & (sub_df['nT_complexity'] <= 6)]
        high_df = sub_df[sub_df['nT_complexity'] >= 7]
        low_mid = sub_df[(sub_df['nT_complexity'] >= 1) & (sub_df['nT_complexity'] <= 6)]
        print(f"{score_col} low (n={low_df.shape[0]}): {low_df[score_col].mean()}")
        print(f"{score_col} mid (n={mid_df.shape[0]}): {mid_df[score_col].mean()}")
        print(f"{score_col} high (n={high_df.shape[0]}): {high_df[score_col].mean()}")
        print(f"{score_col} low_mid (n={low_mid.shape[0]}): {low_mid[score_col].mean()}")

    print('\nSummary statistics by complexity groups (in the paper, we divided into two groups, low(low-mid, <=6) and high (>=7)')
    output_score_for_each_complexity(jaccard_df, 'nT_jaccard')
    output_score_for_each_complexity(precision_df, 'nT_precision')
    output_score_for_each_complexity(recall_df, 'nT_recall')


def compare_aneuploidy_for_nT_clusters(df):
    """
    for clusters without terminal edges, compare the aneuploidy
    for clusters with terminal edges, return 0s
    :param df:
    :return:
    """

    def iterative_nT_aneuploidy(df_row):
        if len(df_row['T_event_ids']):
            return {'TN': 0, 'TP': 0, 'FP': 0, 'FN': 0}
        else:
            file_prefix = data_folder
            file_suffix = '.txt'
            cluster_file_path = file_prefix + df_row['file_name'] + 'cluster_' + str(df_row['cluster']) + file_suffix
            return compare_aneuplodies(cluster_file_path)

    df['nT_aneuploidy_comparison'] = df.apply(lambda row: iterative_nT_aneuploidy(row), axis=1)
    return df


def compare_aneuploidy_for_TnT_clusters(df):
    """
    compare the aneuploidy for all
    :param df:
    :return:
    """

    def iterative_all_aneuploidy(df_row):
        file_prefix = data_folder
        file_suffix = '.txt'
        cluster_file_path = file_prefix + df_row['file_name'] + 'cluster_' + str(df_row['cluster']) + file_suffix
        return compare_aneuplodies(cluster_file_path)

    df['TnT_aneuploidy_comparison'] = df.apply(lambda row: iterative_all_aneuploidy(row), axis=1)
    return df


def summary_statistics_by_event_types(df, status_column='nT_event_status'):
    def rowwise_count_by_event_type(df_row):
        for event_id, event_info in df_row[status_column].items():
            event_type = event_id.split('#')[0]
            if event_type in initial_counts:
                initial_counts[event_type] += event_info['initial count']
                uncaught_counts[event_type] += event_info['count']
            else:
                initial_counts[event_type] = event_info['initial count']
                uncaught_counts[event_type] = event_info['count']

    initial_counts = {}
    uncaught_counts = {}
    df.apply(lambda row: rowwise_count_by_event_type(row), axis=1)

    total_initial = 0
    total_caught = 0
    print(f'\nsummary statistics by event types ({status_column})')
    for event_type, initial_count in initial_counts.items():
        uncaught_count = uncaught_counts[event_type]
        caught_count = initial_count - uncaught_count
        print(f"{event_type}: total={initial_count}, caught={caught_count}, recall={round(caught_count/initial_count, 3)}")
        total_initial += initial_count
        total_caught += caught_count
    print(f"total_initial={total_initial}, total_caught={total_caught}, total_recall={round(total_caught/total_initial, 3)}")


def total_cluster_with_wrong_number_of_chrom(df, aneuploidy_stat_col):
    cluster_with_aneuploidy = [0]
    cluster_without_aneuploidy = [0]
    cluster_with_aneuploidy_wrong = [0]
    cluster_without_aneuploidy_wrong = [0]

    def count_row_has_wrong_number_of_chrom(df_row):
        chrom_status = df_row[aneuploidy_stat_col]
        if chrom_status['TP'] + chrom_status['FN'] > 0:
            cluster_has_aneuploidy = True
            cluster_with_aneuploidy[0] += 1
        else:
            cluster_has_aneuploidy = False
            cluster_without_aneuploidy[0] += 1
        if chrom_status['FP'] + chrom_status['FN'] > 0:
            if cluster_has_aneuploidy:
                cluster_with_aneuploidy_wrong[0] += 1
            else:
                cluster_without_aneuploidy_wrong[0] += 1

    df.apply(lambda row: count_row_has_wrong_number_of_chrom(row), axis=1)
    cluster_with_aneuploidy = cluster_with_aneuploidy[0]
    cluster_without_aneuploidy = cluster_without_aneuploidy[0]
    cluster_with_aneuploidy_wrong = cluster_with_aneuploidy_wrong[0]
    cluster_without_aneuploidy_wrong = cluster_without_aneuploidy_wrong[0]
    print(f'total_cluster_with_wrong_number_of_chrom (n={len(df)})')
    print(f"cluster_with_aneuploidy: {cluster_with_aneuploidy}")
    print(f"cluster_without_aneuploidy: {cluster_without_aneuploidy}")
    print(f"cluster_with_aneuploidy_wrong: {cluster_with_aneuploidy_wrong}")
    print(f"cluster_without_aneuploidy_wrong: {cluster_without_aneuploidy_wrong}")


def compute_row_recall_precision_jaccard(df, score_name, TP, FP, FN):
    df[score_name + '_precision'] = np.where(TP + FP == 0, -1.0, TP / (TP + FP))
    df[score_name + '_recall'] = np.where(TP + FN == 0, -1.0, TP / (TP + FN))
    df[score_name + '_jaccard'] = np.where(TP + FP + FN == 0, -1.0, TP / (TP + FP + FN))
    df[score_name + '_precision'] = df[score_name + '_precision'].round(3)
    df[score_name + '_recall'] = df[score_name + '_recall'].round(3)
    df[score_name + '_jaccard'] = df[score_name + '_jaccard'].round(3)


def iterative_label_missed_SV_edges(df_row):
    # TODO: refactor, this is no longer necessary
    karsim_filename = df_row['file_name']
    karsim_history_edges_filepath = karsim_history_edges_folder + karsim_filename + '.history_sv.txt'
    event_sv_edges = read_history_edges_intermediate_file(karsim_history_edges_filepath)
    missed_sv_edges = df_row['karsim_missed_transition']
    labeled_event_type = []
    # print(event_sv_edges)
    for missed_sv_edge in missed_sv_edges:
        event_found = False
        for entry in event_sv_edges:
            event_type = entry[0]
            event_edges = entry[1]
            if missed_sv_edge in event_edges:
                labeled_event_type.append(event_type)
                event_found = True
                break
        if not event_found:
            labeled_event_type.append('ENF')
    return labeled_event_type


def iterative_get_cn_with_bins(df_row, cn_file_name=get_metadata_file_path('cn_bins_50kbp.txt')):
    cn_bins = read_cn_bin_file(cn_file_name)
    graph = df_row['unmodified_graph']
    karsim_cn, omkar_cn = graph_assign_cn_bin(graph, cn_bins)
    return karsim_cn, omkar_cn


#################CN###############
def read_cn_bin_file(cn_file_name):
    cn_bins = []
    with open(cn_file_name) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            cn_bins.append({'chrom': line[0],
                            'start': int(line[1]),
                            'end': int(line[2])})
    return cn_bins


def filter_cn_bins(chroms, input_cn_bins):
    output_cn_bins = {chrom: [] for chrom in chroms}
    for bin in input_cn_bins:
        if bin['chrom'] in output_cn_bins:
            output_cn_bins[bin['chrom']].append(bin)
    return output_cn_bins


def graph_assign_cn_bin(input_graph: Graph, input_cn_bins):
    ## filter cn_bins to only include Chrs in the current graph
    chroms = input_graph.report_chr_in_graph()
    filtered_cn_bins = filter_cn_bins(chroms, input_cn_bins)
    ## input graphs has disjoint segments
    karsim_cn = {chrom: np.array([0.0 for _ in filtered_cn_bins[chrom]]) for chrom in chroms}
    omkar_cn = {chrom: np.array([0.0 for _ in filtered_cn_bins[chrom]]) for chrom in chroms}
    for start_node, end_nodes in input_graph.karsim_dict.items():
        for end_node in end_nodes:
            edge_type = end_node[2]
            if edge_type == 'segment':
                c_seg = Segment(end_node[0], start_node[1], end_node[1])
                karsim_cn[end_node[0]] += c_seg.assign_cn_bin(filtered_cn_bins[end_node[0]])
    for start_node, end_nodes in input_graph.omkar_dict.items():
        for end_node in end_nodes:
            edge_type = end_node[2]
            if edge_type == 'segment':
                c_seg = Segment(end_node[0], start_node[1], end_node[1])
                omkar_cn[end_node[0]] += c_seg.assign_cn_bin(filtered_cn_bins[end_node[0]])
    return karsim_cn, omkar_cn


def cn_bin_value_triple_state_conversion(input_cn, expected_counts, rounding_allowance=0.05):
    """
    a bin will be called '1'/changed if deviated from WT
    :param rounding_allowance:
    :param expected_counts: WT count for each Chromosome Group
    :param input_cn:
    :return:
    """
    output_triple_state_cn = []
    total_diff = 0
    for chrom, cn_bins in input_cn.items():
        expected_count = expected_counts[chrom]
        for bin_itr in cn_bins:
            diff = bin_itr - expected_count
            total_diff += abs(diff)
            if diff >= rounding_allowance:
                output_triple_state_cn.append(1)
            elif diff <= -1 * rounding_allowance:
                output_triple_state_cn.append(-1)
            else:
                output_triple_state_cn.append(0)
    return np.array(output_triple_state_cn), total_diff


def cn_jaccard_similarity(cn1, cn2):
    intersection = np.sum(np.logical_and(cn1, cn2))
    union = np.sum(np.logical_or(cn1, cn2))
    if union == 0:
        return 1.0  # If both arrays are all zeros, define Jaccard similarity as 1
    else:
        return intersection / union


def average_CN(input_cn_bins):
    output = {}
    for chrom, cn_bins in input_cn_bins.items():
        average_val = int(cn_bins.mean())
        output[chrom] = average_val
    return output


def rowwise_cn_jaccard_sim_by_cluster(df_row, rounding_allowance=0.05):
    # these two dicts need to be identically sized, with each lst in it, also identically sized (pairwise)
    cn_bins1 = df_row['karsim_CN']  # [float]
    cn_bins2 = df_row['omkar_CN']
    # WT expected count is computed using truth genome's average CN
    expected_counts = average_CN(cn_bins1)
    file1_bool_cn, diff1 = cn_bin_value_triple_state_conversion(cn_bins1, expected_counts, rounding_allowance)
    file2_bool_cn, diff2 = cn_bin_value_triple_state_conversion(cn_bins2, expected_counts, rounding_allowance)
    if diff1 <= 1 and diff2 <= 1:
        return 1.0
    else:
        return cn_jaccard_similarity(file1_bool_cn, file2_bool_cn)


##############################################
################ACCURACY BY SV TYPES##########

def count_edges(status_dict):
    total = 0
    for label, info in status_dict.items():
        total += info[0]
    return total

################################################
################ACCURACY BY COMPLEXITY##########

def complexity_by_event_status(input_event_status):
    total_complexity = 0
    for event_id, status in input_event_status.items():
        total_complexity += status['initial count']
    return total_complexity


################################################
################BED FILE PROCESSING#############

def rowwise_match_residual_edge_in_bed_file(df_row, confidence_threshold=0.5):
    def approx_edge_search_bed(input_df, pos1, pos2, ignore_orientation=True, dist=50000):
        if ignore_orientation:
            mask1 = abs(input_df['start'] - pos1) + abs(input_df['end'] - pos2) < dist
            mask2 = abs(input_df['start'] - pos2) + abs(input_df['end'] - pos1) < dist
            return input_df[mask1 | mask2]
    ## read corresponding bed file
    c_file_name = df_row['file_name']

    ## get corresponding residual edge coordinate and event type
    karsim_residuals = df_row['graph'].karsim_edge_label  # dict of {(chr1, pos1, chr2, pos2): [event_id]}

    ## attempt matching each residual edge coordinate against the bed df
    bed_df = df_row['BED']
    bed_df = bed_df[bed_df['confidence'] >= confidence_threshold]
    matches = []
    for e, label in karsim_residuals.items():
        chrom = convert_chrom(e[0][3:])
        filtered_bed_df = bed_df[bed_df['chromosome'] == chrom]
        filtered_bed_df = approx_edge_search_bed(filtered_bed_df, e[1], e[3])
        if not filtered_bed_df.empty:
            matches.append((e[0], e[1], e[2], e[3], label))
    return matches


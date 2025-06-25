import os
import sys
import argparse
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

from read_KarSimulator_output import *
from dependent_clusters_processing import *
from KarUtils import *

parser = argparse.ArgumentParser()
parser.add_argument('--omkar_output', type=str, help='path to DIR of OMKar outputs')
parser.add_argument('--omkar_paths_output', type=str, help="path to DIR containing the OMKar path files")
parser.add_argument('--output_dir', type=str, help="path to this program's output")
parser.add_argument('--karsim_dir', type=str, default='../new_data_files/KarSimulator/', help="path to DIR with corresponding karsim files")
parser.add_argument('--forbidden_region_file', type=str, default=get_metadata_file_path("acrocentric_telo_cen.bed"), help="path to forbidden region file")
args = parser.parse_args()

if not args.omkar_output:
    raise RuntimeError()

omkar_output_folder = args.omkar_output
omkar_path_output_folder = args.omkar_paths_output
output_dir = args.output_dir
karsim_folder = args.karsim_dir
forbidden_region_file = args.forbidden_region_file

os.makedirs(output_dir, exist_ok=True)

println_str = ''

with open(f"{output_dir}/../clustering.log", 'w') as fp_write:
    sys.stdout = fp_write
    for file in os.listdir(omkar_path_output_folder):
        file_name = file.split('.')[0]
        omkar_file_path = omkar_path_output_folder + file
        karsim_file_path = karsim_folder + file_name + '.kt.txt'
        omkar_case_folder = omkar_output_folder + '/' + file.split('/')[-1].replace('.txt', '')
        print(file)
        println_str += str(file) + '\n'

        omkar_log_index_to_breakpoint_dict = omkar_log_get_node(omkar_case_folder)
        omkar_log_diff_edges = omkar_log_get_diff_edges(omkar_case_folder)
        omkar_log_diff_edges_coordinates = translate_indexed_edge_to_coordinates(omkar_log_diff_edges, omkar_log_index_to_breakpoint_dict)
        print(omkar_log_diff_edges_coordinates)
        println_str += str(omkar_log_diff_edges_coordinates) + '\n'

        karsim_segment_to_index_dict, karsim_path_list = read_KarSimulator_output_to_path(karsim_file_path, forbidden_region_file)
        omkar_index_to_segment_dict, omkar_path_list = read_OMKar_output_to_path(omkar_file_path, forbidden_region_file)
        omkar_segment_to_index_dict = reverse_dict(omkar_index_to_segment_dict)
        genome_wide_mutual_breaking(karsim_path_list, omkar_path_list)
        form_dependent_clusters_of_two_path_lists(
            karsim_path_list,
            omkar_path_list,
            karsim_segment_to_index_dict,
            omkar_segment_to_index_dict,
            forbidden_region_file,
            output_dir,
            prefix=file_name,
            omkar_modified_edges=omkar_log_diff_edges_coordinates
        )
        print()
        println_str += '\n'

sys.stdout = sys.__stdout__

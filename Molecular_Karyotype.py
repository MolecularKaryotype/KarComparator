import argparse

from dependent_clusters_processing import *
from KarUtils import *


###########################PARTIAL FILE##################################

def fill_pqter_and_create_seg_obj(input_case, forbidden_region_file='Metadata/acrocentric_telo_cen.bed'):
    input_segments = input_case['segments']
    for key, seg in input_segments.items():
        if seg['start'] == 'pter':
            input_segments[key]['start'] = 0
        else:
            input_segments[key]['start'] = int(input_segments[key]['start'].replace(',', ''))
        if seg['end'] == 'qter':
            input_segments[key]['end'] = get_chr_length_from_forbidden_file(seg['chr_origin'], forbidden_region_file)
        else:
            input_segments[key]['end'] = int(input_segments[key]['end'].replace(',', ''))
    seg_idx_to_obj_dict = {}
    for key, seg in input_segments.items():
        seg_idx_to_obj_dict[key] = Segment(seg['chr_origin'], seg['start'], seg['end'])

    input_case['segment_dict'] = seg_idx_to_obj_dict


def fill_remaining_wt_paths(input_case):
    chr_mentioned = set()
    for path in input_case['paths']:
        chr_mentioned.add(path['chr_name'])
    remaining_wt_paths = []
    for i in range(1, 23):
        if f"Chr{i}" not in chr_mentioned:
            remaining_wt_paths.append(f"Chr{i}")
            remaining_wt_paths.append(f"Chr{i}")
    for sex_chr in input_case['sex']:
        if f"Chr{sex_chr.upper()}" not in chr_mentioned:
            remaining_wt_paths.append(f"Chr{sex_chr.upper()}")

    input_case['remaining_wt'] = remaining_wt_paths


def fill_remaining_segs(input_case, forbidden_region_file='Metadata/acrocentric_telo_cen.bed'):
    chr_to_fill = set(input_case['remaining_wt'])
    new_segs = []
    for chromosome in chr_to_fill:
        chrom_len = get_chr_length_from_forbidden_file(chromosome, forbidden_region_file)
        new_segs.append(Segment(chromosome, 0, chrom_len))
    ## a chrom in 'paths' may all be wt, so it will not have a segment documented
    chrom_in_paths = set()
    for path in input_case['paths']:
        chrom_in_paths.add(path['chr_name'])
    chrom_in_segs = set()
    for seg_idx, seg in input_case['segments'].items():
        chrom_in_segs.add(seg['chr_origin'])
    chrom_not_in_segs = chrom_in_paths - chrom_in_segs  # chrom_in_segs is a subset of chrom_in_paths
    for chromosome in chrom_not_in_segs:
        chrom_len = get_chr_length_from_forbidden_file(chromosome, forbidden_region_file)
        new_segs.append(Segment(chromosome, 0, chrom_len))

    input_case['other_segs'] = new_segs


def rename_segs(input_case):
    all_segments = []
    for seg_idx, seg in input_case['segment_dict'].items():
        all_segments.append(seg)
    all_segments += input_case['other_segs']
    all_segments = sorted(all_segments)
    new_segment_idx_to_obj_dict = {}
    for idx, seg in enumerate(all_segments):
        new_segment_idx_to_obj_dict[idx + 1] = seg
    input_case['full_segment_dict'] = new_segment_idx_to_obj_dict

    ## rename the original paths' segments
    rename_dict = {}
    for original_idx, seg in input_case['segment_dict'].items():
        seg_found = False
        for new_idx, seg_find_itr in input_case['full_segment_dict'].items():
            if seg == seg_find_itr:
                rename_dict[original_idx] = new_idx
                seg_found = True
                break
        if not seg_found:
            raise RuntimeError()
    input_case['rename_dict'] = rename_dict
    for path in input_case['paths']:
        if path['wt']:
            continue
        else:
            path_list = path['path']
        new_path = []
        for path_seg in path_list:
            seg_idx = int(path_seg[:-1])
            seg_sign = path_seg[-1]
            new_seg_idx = rename_dict[seg_idx]
            new_path.append(f"{new_seg_idx}{seg_sign}")
        path['path'] = new_path

    ## fill all wt chrom with spanning segs
    def wt_path(chromosome_name):
        segment_obj_to_idx_dict = reverse_dict(new_segment_idx_to_obj_dict)
        # all_segments is sorted
        spanning_segs = []
        for seg_itr in all_segments:
            if seg_itr.chr_name == chromosome_name:
                spanning_segs.append(seg_itr)
        return_wt_path = [f"{segment_obj_to_idx_dict[seg_obj]}+" for seg_obj in spanning_segs]
        return return_wt_path

    def sort_path_key(input_path):
        val = 0.0
        if input_path['chr_name'] == 'ChrX':
            val += 23
        elif input_path['chr_name'] == 'ChrY':
            val += 24
        else:
            val += int(input_path['chr_name'].replace('Chr', ''))
        if input_path['wt']:
            val += 0.5
        return val

    all_paths = []
    for path in input_case['paths']:
        if path['wt']:
            all_paths.append({'chr_name': path['chr_name'],
                              'path': wt_path(path['chr_name']),
                              'wt': True})
        else:
            all_paths.append(path)
    for path_chrom in input_case['remaining_wt']:
        all_paths.append({'chr_name': path_chrom,
                          'path': wt_path(path_chrom),
                          'wt': True})

    all_paths = sorted(all_paths, key=sort_path_key)
    input_case['all_paths'] = all_paths


def output_molecular_karyotype(input_case, output_dir):
    output_str = "Segment\tNumber\tChromosome\tStart\tEnd\n"
    for seg_idx, seg_obj in input_case['full_segment_dict'].items():
        chr_name = seg_obj.chr_name.replace('Chr', '')
        if chr_name.upper() == 'X':
            chr_name = 23
        elif chr_name.upper() == 'Y':
            chr_name = 24
        output_str += f"Segment\t{seg_idx}\t{chr_name}\t{seg_obj.start}\t{seg_obj.end}\n"
    for path_idx, path in enumerate(input_case['all_paths']):
        output_str += f"Path{path_idx + 1} = {' '.join(path['path'])}\n"
    with open(f"{output_dir}/{input_case['case_name']}.txt", 'w') as fp_write:
        fp_write.write(output_str)


def read_partial_path_list(partial_path_file):
    lines = read_file_into_lines(partial_path_file)
    cases = []
    c_line_idx = 0
    while c_line_idx < len(lines):
        if c_line_idx == 0 and lines[c_line_idx][0] != '>':
            raise ValueError()
        if len(lines[c_line_idx]) == 0:
            c_line_idx += 1
        if lines[c_line_idx][0] == '>':
            case_name = lines[c_line_idx][1:]
            c_line_idx += 1
            if lines[c_line_idx] != '#Sex':
                raise ValueError()
            c_line_idx += 1
            sex = lines[c_line_idx]

            c_line_idx += 1
            if lines[c_line_idx] != '#Segments':
                raise ValueError()
            c_line_idx += 1
            segment_list = {}
            while lines[c_line_idx] != '#Paths':
                c_seg = lines[c_line_idx].split('\t')
                segment_list[int(c_seg[0])] = {'chr_origin': c_seg[1],
                                               'start': c_seg[2],
                                               'end': c_seg[3]}
                c_line_idx += 1

            c_line_idx += 1
            paths = []
            while c_line_idx < len(lines) and len(lines[c_line_idx]) != 0:
                path_info = lines[c_line_idx].split('\t')
                if len(path_info) == 2:
                    path = path_info[1].split(' ')
                    wt = False
                else:
                    path = []
                    wt = True
                paths.append({'chr_name': f"Chr{path_info[0][3:].upper()}",
                              'path': path,
                              'wt': wt})
                c_line_idx += 1
            cases.append({'case_name': case_name,
                          'sex': sex,
                          'segments': segment_list,
                          'paths': paths})
    return cases


def generate_molecular_karyotype_from_partial_path_list(partial_path_file, output_dir):
    cases = read_partial_path_list(partial_path_file)
    for case in cases:
        fill_pqter_and_create_seg_obj(case)
        fill_remaining_wt_paths(case)
        fill_remaining_segs(case)
        rename_segs(case)
        output_molecular_karyotype(case, output_dir)

    return cases


############################################################################
###########################TWO MK CLUSTER###################################

def batch_clustering_for_two_mk_files(input_dir1, input_dir2, output_dir, forbidden_region_file, omkar_debug_dir=None):
    """
    :param input_dir1:
    :param input_dir2:
    :param output_dir:
    :param forbidden_region_file:
    :param omkar_debug_dir: this is corresponding to input_dir2's output
    :return:
    """
    # check if two input folder has one-to-one files
    input_dir1_prefixes = set([filename.split('.')[0] for filename in os.listdir(input_dir1)])
    input_dir2_prefixes = set([filename.split('.')[0] for filename in os.listdir(input_dir2)])
    if input_dir1_prefixes != input_dir2_prefixes:
        raise ValueError("two DIR's files are not one-to-one")

    for filename1 in os.listdir(input_dir1):
        # find filename2
        filename2 = ''
        for filename2_itr in os.listdir(input_dir2):
            if filename1.split('.')[0] == filename2_itr.split('.')[0]:
                filename2 = filename2_itr
                break
        if filename2 == '':
            raise RuntimeError()

        print(f"{filename1}\t{filename2}")
        idx_to_obj_dict1, path_list1 = read_OMKar_output_to_path(f"{input_dir1}/{filename1}", forbidden_region_file)
        idx_to_obj_dict2, path_list2 = read_OMKar_output_to_path(f"{input_dir2}/{filename2}", forbidden_region_file)
        obj_to_idx_dict1 = reverse_dict(idx_to_obj_dict1)
        obj_to_idx_dict2 = reverse_dict(idx_to_obj_dict2)
        genome_wide_mutual_breaking(path_list1, path_list2)
        file_prefix = filename1.split('.')[0]
        if omkar_debug_dir is not None:
            omkar_log_index_to_breakpoint_dict = omkar_log_get_node(omkar_debug_dir)
            omkar_log_diff_edges = omkar_log_get_diff_edges(omkar_debug_dir)
            omkar_log_diff_edges_coordinates = translate_indexed_edge_to_coordinates(omkar_log_diff_edges, omkar_log_index_to_breakpoint_dict)
            print(omkar_log_diff_edges_coordinates)
            form_dependent_clusters_of_two_path_lists(path_list1, path_list2, obj_to_idx_dict1, obj_to_idx_dict2,
                                                      forbidden_region_file, output_dir, omkar_modified_edges=omkar_log_diff_edges_coordinates, prefix=file_prefix)
        else:
            form_dependent_clusters_of_two_path_lists(path_list1, path_list2, obj_to_idx_dict1, obj_to_idx_dict2,
                                                      forbidden_region_file, output_dir, prefix=file_prefix)
        print()


############################################################################


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Molecular Karyotype Utils")
    subparsers = parser.add_subparsers(dest="mode", help="choose a util")

    partial_annotation_parser = subparsers.add_parser("partial_annotation", help="create MK from partial annotation")
    partial_annotation_parser.add_argument("partial_annotation_file", help="file path to the input file")
    partial_annotation_parser.add_argument("output_dir", help="dir path to store the output MK files")

    cluster_mks_parser = subparsers.add_parser("cluster_mks", help="create cluster files between two DIR, both are MKs")
    cluster_mks_parser.add_argument("input_dir1", help="DIR1")
    cluster_mks_parser.add_argument("input_dir2", help="DIR2")
    cluster_mks_parser.add_argument("output_dir", help="output DIR")
    cluster_mks_parser.add_argument("--omkar_debug", type=str, default=None, help="output DIR")
    cluster_mks_parser.add_argument("--forbidden_region_file", type=str, default='Metadata/acrocentric_telo_cen.bed', help="output DIR")

    args = parser.parse_args()
    if args.mode == "partial_annotation":
        os.makedirs(args.output_dir, exist_ok=True)
        generate_molecular_karyotype_from_partial_path_list(args.partial_annotation_file,
                                                            args.output_dir)
    elif args.mode == "cluster_mks":
        os.makedirs(args.output_dir, exist_ok=True)
        batch_clustering_for_two_mk_files(args.input_dir1,
                                          args.input_dir2,
                                          args.output_dir,
                                          args.forbidden_region_file,
                                          args.omkar_debug)

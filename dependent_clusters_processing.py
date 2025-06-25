from KarUtils import *


def genome_wide_mutual_breaking(karsim_path_list, omkar_path_list):
    """
    in a complete bipartite graph fashion perform mutual breaking
    :param karsim_path_list:
    :param omkar_path_list:
    :return:
    """
    for karsim_path in karsim_path_list:
        for omkar_path in omkar_path_list:
            karsim_path.generate_mutual_breakpoints(omkar_path, mutual=True)


def mutually_remove_small_segments(karsim_path_list, omkar_path_list, d=1):
    """
    all segments that are <= d will be removed from both path_lists
    :param karsim_path_list:
    :param omkar_path_list:
    :param d:
    :return:
    """
    # TODO: make this merge instead of direct removal, right now for d=1, we do not cause too big of an impact as we remove
    small_segments = []
    for path in karsim_path_list:
        for segment in path.linear_path.segments:
            if len(segment) <= d:
                if segment not in small_segments:
                    small_segments.append(segment)
    for path in omkar_path_list:
        for segment in path.linear_path.segments:
            if len(segment) <= d:
                if segment not in small_segments:
                    small_segments.append(segment)

    omkar_remove_idxs = []
    karsim_remove_idxs = []
    for path in omkar_path_list:
        path_removal_list = []
        for idx, segment in enumerate(path.linear_path.segments):
            if segment in small_segments:
                path_removal_list.append(idx)
        omkar_remove_idxs.append(path_removal_list)
    for path in karsim_path_list:
        path_removal_list = []
        for idx, segment in enumerate(path.linear_path.segments):
            if segment in small_segments:
                path_removal_list.append(idx)
        karsim_remove_idxs.append(path_removal_list)

    for path_idx, path in enumerate(omkar_path_list):
        current_segments = path.linear_path.segments
        new_segments = [elem for idx, elem in enumerate(current_segments) if idx not in omkar_remove_idxs[path_idx]]
        path.linear_path.segments = new_segments
    for path_idx, path in enumerate(karsim_path_list):
        current_segments = path.linear_path.segments
        new_segments = [elem for idx, elem in enumerate(current_segments) if idx not in karsim_remove_idxs[path_idx]]
        path.linear_path.segments = new_segments


def remove_small_segments(segment_list, d=1):
    remove_idxs = []
    for idx, segment in enumerate(segment_list):
        if len(segment) <= d:
            remove_idxs.append(idx)
    new_segment_list = [elem for idx, elem in enumerate(segment_list) if idx not in remove_idxs]
    return new_segment_list


def extract_cluster_indexed_segments(all_indexed_segments, cluster_chr):
    """
    given a full indexed segment dict, extract only the segments from the list of chr_of_interest
    :param all_indexed_segments:
    :param cluster_chr:
    :return:
    """
    cluster_indexed_segments = []
    for segment in all_indexed_segments:
        if segment.chr_name in cluster_chr:
            cluster_indexed_segments.append(segment)
    return cluster_indexed_segments


def form_dependent_clusters_of_two_path_lists(karsim_path_list,
                                              omkar_path_list,
                                              karsim_segment_to_index_dict,
                                              omkar_segment_to_index_dict,
                                              forbidden_file,
                                              output_dir,
                                              prefix='',
                                              omkar_modified_edges=None):
    ## mutually break segments
    genome_wide_mutual_breaking(karsim_path_list, omkar_path_list)
    ## remove all segment with length <= 1: prevent bug in downstream graph forming
    mutually_remove_small_segments(karsim_path_list, omkar_path_list, d=1)

    ## mutually break the indexing and mark their forbidden region info
    karsim_index_segment_path = Path(Arm(list(karsim_segment_to_index_dict.keys()), 'index_segments'))
    omkar_index_segment_path = Path(Arm(list(omkar_segment_to_index_dict.keys()), 'index_segments'))
    label_path_with_forbidden_regions([karsim_index_segment_path], forbidden_file)
    label_path_with_forbidden_regions([omkar_index_segment_path], forbidden_file)
    karsim_index_segment_path.generate_mutual_breakpoints(omkar_index_segment_path, mutual=True)

    # add forbidden region segments to include the telomeres
    forbidden_region_segment_arm = read_forbidden_regions(forbidden_file)
    forbidden_region_segment_path = Path(forbidden_region_segment_arm)
    karsim_index_segment_path.generate_mutual_breakpoints(forbidden_region_segment_path, mutual=True)
    omkar_index_segment_path.generate_mutual_breakpoints(forbidden_region_segment_path, mutual=True)

    all_indexed_segments = set()
    for segment in karsim_index_segment_path.linear_path.segments:
        all_indexed_segments.add(segment)
    for segment in omkar_index_segment_path.linear_path.segments:
        all_indexed_segments.add(segment)
    for segment in forbidden_region_segment_path.linear_path.segments:
        all_indexed_segments.add(segment)
    all_indexed_segments = sorted(list(all_indexed_segments))

    # remove all segments that have len <= 1
    all_indexed_segments = remove_small_segments(all_indexed_segments, d=1)

    karsim_origins = []
    omkar_origins = []

    for path in karsim_path_list:
        karsim_origins.append(path.get_origins())
    for path in omkar_path_list:
        omkar_origins.append(path.get_origins())

    ## form dependent clusters
    dependent_clusters = form_least_disjoint_supergroups(omkar_origins + karsim_origins)
    clustered_origin_chr = {}
    clustered_karsim_path = {}
    clustered_omkar_path = {}
    for cluster_ind in range(len(dependent_clusters)):
        clustered_origin_chr[cluster_ind] = dependent_clusters[cluster_ind]

        new_karsim_cluster = []
        for path in karsim_path_list:
            path_origin_chr = list(path.get_origins())
            if path_origin_chr[0] in dependent_clusters[cluster_ind]:
                # FIXME: unordered path_origin_chr should equal to dependent_clusters[cluster_ind]
                new_karsim_cluster.append(path)
        clustered_karsim_path[cluster_ind] = new_karsim_cluster

        new_omkar_cluster = []
        for path in omkar_path_list:
            path_origin_chr = list(path.get_origins())
            if path_origin_chr[0] in dependent_clusters[cluster_ind]:
                new_omkar_cluster.append(path)
        clustered_omkar_path[cluster_ind] = new_omkar_cluster

    for cluster_ind in range(len(clustered_origin_chr)):
        current_chr = clustered_origin_chr[cluster_ind]
        current_indexed_segments = extract_cluster_indexed_segments(all_indexed_segments, current_chr)

        output_str = ""
        # cluster info + segment indexing
        output_str += "cluster<{}>, origin_chr<{}>, n_karsim<{}>, n_omkar<{}>\n".format(cluster_ind,
                                                                                        clustered_origin_chr[cluster_ind],
                                                                                        len(clustered_karsim_path[cluster_ind]),
                                                                                        len(clustered_omkar_path[cluster_ind]))
        current_segment_to_index_dict = {}  # also form dict for later parsing
        for segment_ind in range(len(current_indexed_segments)):
            output_str += "{}\t{}\t{}\t{}\t{}\t{}\n".format(segment_ind,
                                                            current_indexed_segments[segment_ind].chr_name,
                                                            current_indexed_segments[segment_ind].start,
                                                            current_indexed_segments[segment_ind].end,
                                                            current_indexed_segments[segment_ind].segment_type,
                                                            len(current_indexed_segments[segment_ind]))
            current_segment_to_index_dict[current_indexed_segments[segment_ind]] = str(segment_ind)
        output_str += "---\n"



        # This snippet prevents path still having unbroken seg relative to the dict
        def break_path_by_dictionary(dictionary_path: Path, input_path: Path):
            input_path.generate_mutual_breakpoints(dictionary_path, mutual=False)

        c_segment_to_index_dict_seg_list = []
        for seg in current_segment_to_index_dict:
            c_segment_to_index_dict_seg_list.append(seg)
        c_segment_to_index_dict_path = Path(Arm(c_segment_to_index_dict_seg_list, 'dict_segs'))

        # karsim path
        for path in clustered_karsim_path[cluster_ind]:
            break_path_by_dictionary(c_segment_to_index_dict_path, path)
            output_str += path.tostring_path_by_index(current_segment_to_index_dict) + '\n'
        output_str += "---\n"

        # omkar path
        for path in clustered_omkar_path[cluster_ind]:
            break_path_by_dictionary(c_segment_to_index_dict_path, path)
            output_str += path.tostring_path_by_index(current_segment_to_index_dict) + '\n'
        output_str += "---\n"

        # include omkar introduced edges
        if omkar_modified_edges is not None:
            for edge, multiplicity in omkar_modified_edges.items():
                node1_chr = edge[0]
                node1_coor = edge[1]
                node2_chr = edge[2]
                node2_coor = edge[3]
                edge_type = edge[4]

                if (node1_chr in current_chr) != (node2_chr in current_chr):
                    raise RuntimeError('only one breakpoint of the edge in current dependent component')
                if node1_chr in current_chr and node2_chr in current_chr:
                    # edge is in current dependent component, index it to node endpoints
                    node1_idx = None
                    node2_idx = None
                    for segment, segment_idx in current_segment_to_index_dict.items():
                        if node1_chr == segment.chr_name:
                            if node1_coor == segment.start:
                                node1_idx = str(segment_idx) + 's'
                            elif node1_coor == segment.end:
                                node1_idx = str(segment_idx) + 't'
                        if node2_chr == segment.chr_name:
                            if node2_coor == segment.start:
                                node2_idx = str(segment_idx) + 's'
                            elif node2_coor == segment.end:
                                node2_idx = str(segment_idx) + 't'
                    if node1_idx is None or node2_idx is None:
                        raise RuntimeError('edge labeling for OMKar diff edge unsuccessful: ' + str(edge))
                    if node1_chr == node2_chr:
                        edge_distance = abs(node1_coor - node2_coor)
                    else:
                        edge_distance = -1
                    output_str += '({}, {}, {}, {})\t{}\n'.format(node1_idx, node2_idx, multiplicity, edge_type, edge_distance)
            output_str += "---\n"

        with open("{}/{}cluster_{}.txt".format(output_dir, prefix, str(cluster_ind)), 'w') as fp_write:
            fp_write.write(output_str)


def form_least_disjoint_supergroups(all_origins):
    """
    Given a list of values, form supergroups that are the smallest size but disjoint in the initial grouping
    :param all_origins:
    :return:
    """

    class Bin:
        def __init__(self, chr_components: set):
            self.chr_components = set(chr_components)

        def __len__(self):
            return len(self.chr_components)

        def __str__(self):
            return str(self.chr_components)

        def includes_chr(self, input_chr):
            return input_chr in self.chr_components

        def chr_add(self, input_chr):
            self.chr_components.add(input_chr)

        def merge_bins(self, other):
            new_components = self.chr_components
            for new_chr in other.chr_components:
                new_components.add(new_chr)
            new_bin = Bin(new_components)
            return new_bin

        def get_sorted_chr_list(self):
            chr_list = list(self.chr_components)

            def custom_sort_key(chromosome):
                numerical_value = chromosome[3:]
                if numerical_value.upper() == 'X':
                    numerical_value = 23
                elif numerical_value.upper() == 'Y':
                    numerical_value = 24
                else:
                    numerical_value = int(numerical_value)
                return numerical_value

            chr_list = sorted(chr_list, key=custom_sort_key)
            return chr_list

    bin_list = []
    for origin in all_origins:
        bin_to_construct = Bin(origin)
        for origin_chr in origin:
            # if an intersection found in a bin, add all members of the group into that bin
            for bin_itr in bin_list:
                if bin_itr.includes_chr(origin_chr):
                    bin_to_construct = bin_to_construct.merge_bins(bin_itr)

        # search rest of the chr in the merged bin, so we can merge their corresponding bins also
        bin_to_remove = []
        for bin_ind in range(len(bin_list)):
            current_bin_itr = bin_list[bin_ind]
            if current_bin_itr.chr_components.intersection(bin_to_construct.chr_components):
                bin_to_remove.append(bin_ind)
                bin_to_construct = bin_to_construct.merge_bins(current_bin_itr)

        if len(bin_to_construct) == 0:
            # none of the origin chrs is in a bin yet
            bin_to_construct = Bin(origin)
        else:
            new_bin_list = [value for index, value in enumerate(bin_list) if index not in bin_to_remove]
            bin_list = new_bin_list

        bin_list.append(bin_to_construct)

    output = []
    for bin_itr in bin_list:
        output.append(bin_itr.get_sorted_chr_list())
    return output


def omkar_log_get_node(folder_path):
    """
    Get node_to_breakpoint dictionary
    :param folder_path: omkar output's case folder
    :return:
    """
    node_file = folder_path + "/" + folder_path.split('/')[-1] + '.preILP_nodes.txt'
    index_to_breakpoint_dict = {}
    with open(node_file) as fp_read:
        for line in fp_read:
            info = line.replace('\n', '').split('\t')
            chr_origin = info[1]
            if chr_origin == '23':
                chr_origin = "X"
            elif chr_origin == '24':
                chr_origin = "Y"
            chr_origin = 'Chr' + chr_origin
            pos = round_half_up(float(info[2]))
            index_to_breakpoint_dict[info[0]] = (chr_origin, pos)
    return index_to_breakpoint_dict


def omkar_log_get_diff_edges(folder_path):
    """
    Get all the dummy edge AND SEG edge that changed value in this whole case
    :param folder_path: omkar output's case folder
    :return:
    """
    pre_ILP_edges = {}
    post_ILP_edges = {}

    pre_ILP_edge_file = folder_path + "/" + folder_path.split('/')[-1] + '.preILP_edges.txt'
    with open(pre_ILP_edge_file) as fp_read:
        for line in fp_read:
            info = line.replace('\n', '').replace('(', '').replace(')', '').split(', ')
            info[3] = info[3].replace("'", '')
            pre_ILP_edges[(info[0], info[1], info[3])] = int(info[2])

    for dummy_file in os.listdir(folder_path + '/all_edges_with_dummies'):
        with open(folder_path + '/all_edges_with_dummies/' + dummy_file) as fp_read:
            for line in fp_read:
                info = line.replace('\n', '').replace('(', '').replace(')', '').split(', ')
                info[3] = info[3].replace("'", '')
                post_ILP_edges[(info[0], info[1], info[3])] = int(info[2])

    ## compute diff: only compare SEG edge and append Dummy edge
    diff_edges = {}
    for edge, post_ILP_count in post_ILP_edges.items():
        if edge[2] == 'S':
            if edge not in pre_ILP_edges:
                diff_edges[edge] = post_ILP_count
            else:
                pre_ILP_count = pre_ILP_edges[edge]
                if post_ILP_count - pre_ILP_count != 0:
                    diff_edges[edge] = post_ILP_count - pre_ILP_count
                pre_ILP_edges.pop(edge)
        elif edge[2] == 'D':
            diff_edges[edge] = post_ILP_count
    for remain_edge, pre_ILP_count in pre_ILP_edges.items():
        if remain_edge[2] == 'S':
            if remain_edge in diff_edges:
                raise RuntimeError('double add, pop not working properly')
            if pre_ILP_count != 0:
                diff_edges[remain_edge] = -1 * pre_ILP_count

    return diff_edges


def translate_indexed_edge_to_coordinates(edge_list: {(str, str, str): int},
                                          index_to_breakpoint_dict: {str: (str, int)}) -> {(str, int, str, int, str): int}:
    """
    :param edge_list: {(node1_idx, node2_idx, edge_type): multiplicity}
    :param index_to_breakpoint_dict: {node_idx: (chr_origin, breakpoint_coordinate)}
    :return: {(chr1, coor1, chr2, coor2, edge_type): multiplicity}
    """
    return_dict = {}
    for edge, multiplicity in edge_list.items():
        edge_type = edge[2]
        node1 = index_to_breakpoint_dict[edge[0]]
        node2 = index_to_breakpoint_dict[edge[1]]
        return_dict[(node1[0], node1[1], node2[0], node2[1], edge_type)] = multiplicity
    return return_dict
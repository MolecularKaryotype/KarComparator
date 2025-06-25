import os

import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
import argparse
import ast

from KarUtils import read_forbidden_regions


forbidden_region_file = '../KarUtils/Metadata/acrocentric_telo_cen.bed'
uniform_dist = 1/6


class Vertex:
    name: str
    origin_chr: str
    pos: int

    def __init__(self, name, origin_chr, pos):
        self.name = name
        self.origin_chr = origin_chr
        self.pos = pos

    def __str__(self):
        return '<name: {}, chr: {}, pos: {}>'.format(self.name, self.origin_chr, self.pos)


def get_vertices_pre_ILP(log_file):
    V = {}
    with open(log_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            pos = int(line[2].split('.')[0])
            V[line[0]] = Vertex(int(line[0]), line[1], pos)
    return V


def rename_and_filter_vertices(current_V, new_name_dict, chr_of_interest):
    new_V = {}
    rename_dict = {}
    for old_name, vertex in current_V.items():
        old_chr = vertex.origin_chr
        old_pos = vertex.pos
        if old_chr not in chr_of_interest:
            continue
        if old_chr == '23':
            new_chr = 'ChrX'
        elif old_chr == '24':
            new_chr = 'ChrY'
        else:
            new_chr = 'Chr' + str(old_chr)
        new_name = new_name_dict[(new_chr, old_pos)]
        new_vertex = Vertex(new_name, new_chr, old_pos)
        new_V[new_name] = new_vertex
        rename_dict[old_name] = new_name
    return new_V, rename_dict


def get_edges_pre_ILP(log_file):
    E = []
    with open(log_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').replace('(', '').replace(')', '').replace("'", '').split(', ')
            E.append(tuple(line))
    return E


def get_edges_post_ILP(log_file):
    E = []
    with open(log_file) as fp_read:
        fp_read.readline()  # skip first line
        for line in fp_read:
            line = line.replace('\n', '').replace('(', '').replace(')', '').replace("'", '').split(', ')
            E.append(tuple(line))
    return E


def get_dummy_edges(log_file):
    E = []
    with open(log_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').replace('(', '').replace(')', '').replace("'", '').split(', ')
            if line[3] == 'D':
                # only include the dummy edges
                E.append(tuple(line))
    return E


def consolidate_edge_with_dummy_edge(E, dummy_E):
    for dummy_E_itr in dummy_E:
        val0 = dummy_E_itr[0]
        val1 = dummy_E_itr[1]

        same_edge_found = False
        for E_itr in E:
            if E_itr[0] == val0 and E_itr[1] == val1:
                same_edge_found = True
                break

        if not same_edge_found:
            E.append(dummy_E_itr)
        else:
            # test if the reverse of the dummy edge is also in E
            for E_itr in E:
                if E_itr[0] == val1 and E_itr[1] == val0:
                    raise RuntimeError('reverse of dummy edge also in file, preventing overwrite')
            reversed_dummy = tuple([dummy_E_itr[1], dummy_E_itr[0], dummy_E_itr[2], dummy_E_itr[3]])
            E.append(reversed_dummy)
    return E


def get_cluster_of_interest(current_V, chrs_of_interest, cluster_metadata):
    filtered_V_names = set()
    for v_name, vertex in current_V.items():
        v_chr = vertex.origin_chr
        if v_chr not in chrs_of_interest:
            continue
        else:
            filtered_V_names.add(v_name)
    cluster_number = set()
    for cluster_idx, cluster_nodes in cluster_metadata.items():
        if filtered_V_names.intersection(cluster_nodes):
            cluster_number.add(cluster_idx)
    return list(cluster_number)


def filter_nodes_and_edges(current_E, current_V, chrs_of_interest):
    """
    given the chrs_of_interest, get all nodes with these origin_chr;
    then, for all edges invovled with these nodes, add both nodes from the edge and include the edge
    :param current_E:
    :param current_V:
    :param chrs_of_interest:
    :return:
    """
    filtered_V_names = set()
    filtered_V = {}
    filtered_E = []
    for v_name, vertex in current_V.items():
        v_chr = vertex.origin_chr
        if v_chr not in chrs_of_interest:
            continue
        else:
            filtered_V_names.add(v_name)

    for edge in current_E:
        v1 = edge[0]
        v2 = edge[1]
        if v1 in filtered_V_names or v2 in filtered_V_names:
            filtered_V_names.add(v1)
            filtered_V_names.add(v2)
            filtered_E.append(edge)

    for v_name in filtered_V_names:
        filtered_V[v_name] = current_V[v_name]

    return filtered_E, filtered_V


def rename_and_filter_edges(current_E, rename_dict):
    vertex_of_interest = list(rename_dict.keys())
    new_E = []
    for old_edge in current_E:
        node1 = old_edge[0]
        node2 = old_edge[1]
        if node1 not in vertex_of_interest or node2 not in vertex_of_interest:
            continue
        else:
            new_edge = (rename_dict[node1], rename_dict[node2], old_edge[2], old_edge[3])
            new_E.append(new_edge)
    return new_E


def filter_nodes(chrs_of_interest, V):
    key_to_pop = []
    for key in V:
        if V[key].origin_chr not in chrs_of_interest:
            key_to_pop.append(key)
    for key in key_to_pop:
        V.pop(key)
    return V


def iterative_add_edge(edges_list: [(str, str, str, str)], graph):
    for edge in edges_list:
        node1 = edge[0]
        node2 = edge[1]
        multiplicity = int(edge[2])
        edge_type = edge[3]

        edge_color = None
        if edge_type == 'S':
            edge_color = 'black'
        elif edge_type == 'R':
            edge_color = 'blue'
        elif edge_type == 'SV':
            edge_color = 'red'
        elif edge_type == 'D':
            edge_color = 'orange'

        graph.add_edge(node1, node2, color=edge_color, weight=multiplicity)


def generate_uniform_linear_coordinates(n, fixed_y=0.5, fixed_distance=uniform_dist):
    coordinates = []
    for i in range(n):
        x = round(0 + i * fixed_distance, 4)
        y = fixed_y
        coordinates.append((x, y))

    return coordinates


def label_centromere_nodes(V, forbidden_file):
    """
    archived, should label ref edge as CEN edge,instead of nodes
    :param V:
    :param forbidden_file:
    :return:
    """
    forbidden_region_segments = read_forbidden_regions(forbidden_file).segments
    centromere_segments = []
    for segment in forbidden_region_segments:
        if 'centromere' in segment.segment_type:
            centromere_segments.append(segment)

    centromere_V = []
    for node_id, node in V.items():
        current_chr = 'Chr' + node.origin_chr
        current_chr_centromere_segment = None
        for segment in centromere_segments:
            if segment.chr_name == current_chr:
                current_chr_centromere_segment = segment
                break
        if current_chr_centromere_segment.start <= node.pos <= current_chr_centromere_segment.end:
            centromere_V.append(node_id)
    print('x')


def make_omkar_graph(filtered_V, filtered_E, pic_output_path):
    graph_width = len(filtered_V)
    fixed_height = 6

    G = nx.DiGraph()
    G.add_nodes_from(filtered_V)
    iterative_add_edge(filtered_E, G)

    V_pos = {}
    V_names = list(filtered_V.keys())
    cor = generate_uniform_linear_coordinates(len(filtered_V))
    ordered_V_names = sorted(V_names, key=int)
    for node_ind in range(len(filtered_V)):
        V_pos[ordered_V_names[node_ind]] = cor[node_ind]

    E_colors = nx.get_edge_attributes(G, 'color')
    E_weights = nx.get_edge_attributes(G, 'weight')

    plt.figure(figsize=(graph_width, fixed_height))
    plot = ng.InteractiveGraph(G,
                               node_layout=V_pos,
                               node_labels=True,
                               edge_color=E_colors,
                               edge_layout='arc',
                               edge_labels=E_weights,
                               arrows=True,
                               node_size=6,
                               node_label_offset=0.001,
                               node_label_font_dict=dict(size=10),
                               edge_label_fontdict=dict(size=9),
                               scale=(graph_width, 1))
    plt.savefig(pic_output_path)
    plt.close()


def graph_pre_ILP(chrs_of_interest, header_name, output_dir, V_rename_dict):
    node_file = dir_name + header_name + '/' + header_name + '.preILP_nodes.txt'
    edge_file = dir_name + header_name + '/' + header_name + '.preILP_edges.txt'

    V = get_vertices_pre_ILP(node_file)
    E = get_edges_pre_ILP(edge_file)

    if len(V_rename_dict) != 0:
        # TODO: filtering logic can be done better, see the else case logic
        filtered_V, E_rename_dict = rename_and_filter_vertices(V, V_rename_dict, chrs_of_interest)
        filtered_E = rename_and_filter_edges(E, E_rename_dict)
    else:
        filtered_E, filtered_V = filter_nodes_and_edges(E, V, chrs_of_interest)
    print(filtered_V.keys())

    pic_output_path = output_dir + header_name + str(chrs_of_interest) + '.preILP.png'
    make_omkar_graph(filtered_V, filtered_E, pic_output_path)


def graph_post_ILP(chrs_of_interest, header_name, output_dir, V_rename_dict):
    node_file = dir_name + header_name + '/' + header_name + '.preILP_nodes.txt'
    metadata_file = dir_name + header_name + '/postILP_components/' + header_name + '.postILP.metadata.txt'

    V = get_vertices_pre_ILP(node_file)

    ## find the right cluster file
    cluster_metadata = {}
    with open(metadata_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            cluster_number = line[0]
            cluster_nodes = [str(i) for i in eval(line[1])]
            cluster_metadata[cluster_number] = cluster_nodes
    # assumes that only one cluster is in play
    cluster_to_graph = get_cluster_of_interest(V, chrs_of_interest, cluster_metadata)
    if len(cluster_to_graph) != 1:
        print(cluster_to_graph)
        raise RuntimeError()
    else:
        cluster_to_graph = cluster_to_graph[0]

    edge_file = dir_name + header_name + '/postILP_components/' + header_name + '.postILP_component_' + cluster_to_graph + '.txt'
    E = get_edges_post_ILP(edge_file)

    if len(V_rename_dict) != 0:
        filtered_V, E_rename_dict = rename_and_filter_vertices(V, V_rename_dict, chrs_of_interest)
        filtered_E = rename_and_filter_edges(E, E_rename_dict)
    else:
        filtered_E, filtered_V = filter_nodes_and_edges(E, V, chrs_of_interest)

    pic_output_path = output_dir + header_name + str(chrs_of_interest) + '.postILP.png'
    make_omkar_graph(filtered_V, filtered_E, pic_output_path)


def graph_post_ILP_with_dummies(chrs_of_interest, header_name, output_dir, V_rename_dict):
    node_file = dir_name + header_name + '/' + header_name + '.preILP_nodes.txt'
    metadata_file = dir_name + header_name + '/postILP_components/' + header_name + '.postILP.metadata.txt'

    V = get_vertices_pre_ILP(node_file)

    ## find the right cluster file
    cluster_metadata = {}
    with open(metadata_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            cluster_number = line[0]
            cluster_nodes = [str(i) for i in eval(line[1])]
            cluster_metadata[cluster_number] = cluster_nodes
    # assumes that only one cluster is in play
    cluster_to_graph = get_cluster_of_interest(V, chrs_of_interest, cluster_metadata)
    if len(cluster_to_graph) != 1:
        raise RuntimeError()
    else:
        cluster_to_graph = cluster_to_graph[0]

    edge_file = dir_name + header_name + '/postILP_components/' + header_name + '.postILP_component_' + cluster_to_graph + '.txt'
    dummy_file = dir_name + header_name + '/all_edges_with_dummies/' + header_name + '.with_dummies_component_' + cluster_to_graph + '.txt'
    E = get_edges_post_ILP(edge_file)
    dummy_E = get_dummy_edges(dummy_file)

    ## prevent dummy_edge overwriting regular edge by having duplicate (multi-graph not allowed)
    E = consolidate_edge_with_dummy_edge(E, dummy_E)

    if len(V_rename_dict) != 0:
        filtered_V, E_rename_dict = rename_and_filter_vertices(V, V_rename_dict, chrs_of_interest)
        filtered_E = rename_and_filter_edges(E, E_rename_dict)
    else:
        filtered_E, filtered_V = filter_nodes_and_edges(E, V, chrs_of_interest)

    pic_output_path = output_dir + header_name + str(chrs_of_interest) + '.with_dummy.png'
    make_omkar_graph(filtered_V, filtered_E, pic_output_path)


if __name__ == "__main__":
    def parse_list(s):
        x = ast.literal_eval(s)
        new_list = []
        for item in x:
            new_list.append(str(item))
        return new_list

    parser = argparse.ArgumentParser(description='plotting omkar preILP and postILP graphs')
    parser.add_argument('data_dir', type=str)
    parser.add_argument('case_name', type=str)
    parser.add_argument('output_dir', type=str)
    parser.add_argument('chr_of_int', type=parse_list)
    parser.add_argument('--rename_dict', default='{}', type=str)

    args = parser.parse_args()
    dir_name = args.data_dir
    chr_of_int = args.chr_of_int
    current_output_dir = args.output_dir
    header = args.case_name
    rename_dict = ast.literal_eval(args.rename_dict)

    os.makedirs(current_output_dir, exist_ok=True)

    graph_pre_ILP(chr_of_int, header, current_output_dir, rename_dict)
    graph_post_ILP(chr_of_int, header, current_output_dir, rename_dict)
    graph_post_ILP_with_dummies(chr_of_int, header, current_output_dir, rename_dict)

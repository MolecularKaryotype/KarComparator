import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
import math
from collections import Counter
from scipy.optimize import linear_sum_assignment

from read_cluster_file import *
from read_KarSimulator_output import *


i_matching_allowance = 200000
i_approximate_allowance = 50000


def is_self_edge(chr1, pos1, chr2, pos2):
    if (chr1, pos1) == (chr2, pos2):
        return True
    return False

def intra_transition_edge_distance(chr1, pos1, chr2, pos2):
    """
    Distances between two nodes
    """
    if chr1 != chr2:
        return -1
    else:
        return abs(pos2 - pos1)


def inter_transition_edge_distance(edge1_chr1, edge1_pos1, edge1_chr2, edge1_pos2,
                                   edge2_chr1, edge2_pos1, edge2_chr2, edge2_pos2):
    """
    Distances between two edges
    """
    if edge1_chr1 != edge2_chr1 or edge1_chr2 != edge2_chr2:
        return -1
    # prevent double-count for the dupinv's self edge (we calc distance once)
    if (edge1_chr1, edge1_pos1) == (edge1_chr2, edge1_pos2) and (edge2_chr1, edge2_pos1) == (edge2_chr2, edge2_pos2):
        return abs(edge1_pos1 - edge2_pos1)
    else:
        return abs(edge1_pos1 - edge2_pos1) + abs(edge1_pos2 - edge2_pos2)


def custom_sort_node(node_tuple):
    chr_info = node_tuple[0]
    pos_info = int(node_tuple[1])
    chr_index = chr_info[3:]
    if chr_index == "X":
        chr_index = 23
    elif chr_index == "Y":
        chr_index = 24
    else:
        chr_index = int(chr_index)

    return chr_index, pos_info


def pop_edge(chr1, pos1, chr2, pos2, edge_type, target_dict):
    node2_list = target_dict[(chr1, pos1)]
    if (chr2, pos2, edge_type) not in node2_list:
        raise ValueError('edge does not exist')
    else:
        node2_list.remove((chr2, pos2, edge_type))


class Graph:
    node_name: {(str, int): str}  # chr, position: node name
    source_sink_nodes: {(str, int): str}  # chr, pos: name; used for computing distance
    segment_edge_type: {(str, int, str, int): str}  # documents for comparison forbidden regions
    karsim_dict: {(str, int): [(str, int, str)]}  # chr1, pos1: chr2, pos2, edge_type
    karsim_edge_label: {(str, int, str, int): [str]}  # label of what SV introduced this edge
    omkar_dict: {(str, int): [(str, int, str)]}  # chr1, pos1: chr2, pos2, edge_type
    approximated_cnv: int
    karsim_n_transition_approximated: int
    omkar_n_transition_approximated: int
    # TODO: refactor this to omkar_edge_label
    edges_of_interest: [(str, str, int, str, int)]  # start_node, end_node, multiplicity, edge_type, edge_distance
    events: {str: int}  # event_name: multiplicity

    def __init__(self):
        self.node_name = {}
        self.node_name_reverse_dict = {}
        self.segment_edge_type = {}
        self.karsim_dict = {}
        self.karsim_edge_label = {}
        self.omkar_dict = {}
        self.approximated_cnv = 0
        self.karsim_n_transition_approximated = 0
        self.omkar_n_transition_approximated = 0
        self.edges_of_interest = []
        self.events = {}
        self.translated_karsim_dict = {}
        self.translated_omkar_dict = {}
        self.source_sink_nodes = {}

    def generate_node_name_reverse_dict(self):
        self.node_name_reverse_dict = reverse_dict(self.node_name)

    def add_sv_label_to_karsim_edges(self, event_sv_edges):
        """
        label each Karsim sv edges with event labels
        :param event_sv_edges: a DICT containing event sv labels
        :return:
        """
        labeled_edges = {}
        for start_node, end_nodes in self.karsim_dict.items():
            for end_node in end_nodes:
                if start_node in self.node_name:
                    start_node_orientation = self.node_name[start_node][-1]
                else:
                    start_node_orientation = self.source_sink_nodes[start_node][0]
                if end_node[:2] in self.node_name:
                    end_node_orientation = self.node_name[end_node[:2]][-1]
                else:
                    end_node_orientation = self.source_sink_nodes[end_node[:2]][0]
                edge_orientation = (start_node_orientation, end_node_orientation)
                if end_node[2] == 'segment':
                    continue
                elif end_node[0] == start_node[0] and abs(end_node[1] - start_node[1]) <= 5 \
                        and (edge_orientation == ('t', 's') or edge_orientation == ('s', 't')):
                    # skip reference edges
                    continue
                c_edge = (*start_node, *end_node[:2])
                c_edge_event_type = ''
                for entry in event_sv_edges:
                    entry_found = -1
                    for edge_idx, event_edge in enumerate(entry['edges']):
                        if edges_are_similar(c_edge, event_edge):
                            c_edge_event_type = entry['label']
                            entry_found = edge_idx
                            break
                    if entry_found != -1:
                        if entry['label'] in labeled_edges:
                            labeled_edges[entry['label']].append(entry['edges'].pop(entry_found))
                        else:
                            labeled_edges[entry['label']] = [entry['edges'].pop(entry_found)]
                        break
                if c_edge_event_type == '':
                    c_edge_event_type = 'ENF'
                if c_edge not in self.karsim_edge_label:
                    self.karsim_edge_label[c_edge] = [c_edge_event_type]
                else:
                    self.karsim_edge_label[c_edge].append(c_edge_event_type)
        return labeled_edges

    def tally_karsim_edge_event_types(self):
        """
        this is for sanity checking whether all SV edges are labeled with an event type
        :param current_tally_dict: dict to be appended to
        :return: dict {event_type: multiplicity}
        """
        current_tally_dict = {}
        for edge, sv_event_types in self.karsim_edge_label.items():
            # see if edge still on the graph
            for start_node, end_nodes in self.karsim_dict.items():
                for end_node in end_nodes:
                    c_edge = (*start_node, *end_node[:2])
                    if edge == c_edge:
                        for sv_event_type in sv_event_types:
                            if sv_event_type in current_tally_dict:
                                current_tally_dict[sv_event_type] += 1
                            else:
                                current_tally_dict[sv_event_type] = 1
        return current_tally_dict

    def add_edge_to_dict(self, chr1, pos1, chr2, pos2, edge_type, target_dict: str):
        if target_dict == 'omkar':
            target_dict = self.omkar_dict
        elif target_dict == 'karsim':
            target_dict = self.karsim_dict
        else:
            raise ValueError()

        if (chr1, pos1) in target_dict:
            target_dict[(chr1, pos1)].append((chr2, pos2, edge_type))
        else:
            target_dict[(chr1, pos1)] = [(chr2, pos2, edge_type)]

    def add_segment_edge(self, input_segment: Segment, target_graph: str):
        """
        :param input_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        # add node
        if input_segment.direction():
            self.node_name[(input_segment.chr_name, input_segment.start)] = input_segment.kt_index[:-1] + "s"
            self.node_name[(input_segment.chr_name, input_segment.end)] = input_segment.kt_index[:-1] + "t"
        else:
            self.node_name[(input_segment.chr_name, input_segment.end)] = input_segment.kt_index[:-1] + "s"
            self.node_name[(input_segment.chr_name, input_segment.start)] = input_segment.kt_index[:-1] + "t"

        # document edge type
        self.segment_edge_type[(input_segment.chr_name, input_segment.start,
                                input_segment.chr_name, input_segment.end)] = input_segment.segment_type

        # add edge to dict
        self.add_edge_to_dict(input_segment.chr_name, input_segment.start,
                              input_segment.chr_name, input_segment.end,
                              'segment', target_graph)

    def add_transition_edge(self, start_segment: Segment, end_segment: Segment, target_graph: str):
        """
        transition edge from the end of the start_segment to the start of the end_segment
        :param start_segment:
        :param end_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        # add edge to dict
        self.add_edge_to_dict(start_segment.chr_name, start_segment.end,
                              end_segment.chr_name, end_segment.start,
                              'transition', target_graph)

    def add_source_transition_edge(self, source_pos, first_segment, target_graph: str):
        # source node always have the same chr_name as the first seg
        self.add_edge_to_dict(first_segment.chr_name, source_pos,
                              first_segment.chr_name, first_segment.start,
                              'transition_source', target_graph)

    def add_sink_transition_edge(self, end_pos, last_segment, target_graph: str):
        # sink node always have the same chr_name as the last seg
        self.add_edge_to_dict(last_segment.chr_name, last_segment.end,
                              last_segment.chr_name, end_pos,
                              'transition_sink', target_graph)

    def prune_same_edges(self):
        dict1 = self.karsim_dict
        dict2 = self.omkar_dict
        matched_labels = []

        # Find the common keys between the dictionaries
        common_keys = dict1.keys() & dict2.keys()

        # Iterate over common keys and remove instances
        for node1 in common_keys:
            occurrences_in_dict1 = Counter(dict1[node1])
            occurrences_in_dict2 = Counter(dict2[node1])

            for node2 in occurrences_in_dict1:
                # Determine the minimum number of occurrences between the two dicts
                if node2 in occurrences_in_dict2:
                    # this also means the edge types are the same
                    min_occurrences = min(occurrences_in_dict1[node2], occurrences_in_dict2[node2])
                    for _ in range(min_occurrences):
                        dict1[node1].remove(node2)
                        dict2[node1].remove(node2)
                        if node2[2] != 'segment':
                            unique_positions = {node1[1], node2[1]}
                            if (node1[0], node1[1], node2[0], node2[1]) in self.karsim_edge_label:
                                # 0s for distances as we have perfect matching here
                                matched_labels.append({'event_id': self.karsim_edge_label[(node1[0], node1[1], node2[0], node2[1])][0],
                                                       'distances': [0 for _ in unique_positions],
                                                       'karsim_edge': (node1[0], node1[1], node2[0], node2[1]),
                                                       'omkar_edge': (node1[0], node1[1], node2[0], node2[1])})
                                if len(self.karsim_edge_label[(node1[0], node1[1], node2[0], node2[1])]) == 1:
                                    self.karsim_edge_label.pop((node1[0], node1[1], node2[0], node2[1]))
                                else:
                                    self.karsim_edge_label[(node1[0], node1[1], node2[0], node2[1])].pop(0)
                            else:
                                # 0s for distances as we have perfect matching here
                                matched_labels.append({'event_id': 'ENF',
                                                       'distances': [0 for _ in unique_positions],
                                                       'karsim_edge': (node1[0], node1[1], node2[0], node2[1]),
                                                       'omkar_edge': (node1[0], node1[1], node2[0], node2[1])})

        # Remove keys that have become empty
        self.karsim_dict = {k: v for k, v in dict1.items() if v}
        self.omkar_dict = {k: v for k, v in dict2.items() if v}

        return matched_labels

    def gather_edges(self, target_graph: str) -> ({(str, str): int}, {(str, str): int}):
        """
        get all segment and transition edges
        :param target_graph: 'karsim' or 'omkar'
        :return:
        """
        if target_graph == 'karsim':
            target_dict = self.karsim_dict
        elif target_graph == 'omkar':
            target_dict = self.omkar_dict
        else:
            raise ValueError

        E_segment = {}
        E_transition = {}
        for node1, value in target_dict.items():
            for node2 in value:
                if (node1[0], node1[1]) in self.node_name:
                    node1_name = self.node_name[(node1[0], node1[1])]
                else:
                    node1_name = self.source_sink_nodes[(node1[0], node1[1])]
                if (node2[0], node2[1]) in self.node_name:
                    node2_name = self.node_name[(node2[0], node2[1])]
                else:
                    node2_name = self.source_sink_nodes[(node2[0], node2[1])]
                edge_type = node2[2]

                if edge_type == 'segment':
                    if (node1_name, node2_name) in E_segment:
                        E_segment[(node1_name, node2_name)] += 1
                    else:
                        E_segment[(node1_name, node2_name)] = 1
                elif edge_type.startswith('transition'):
                    if (node1_name, node2_name) in E_transition:
                        E_transition[(node1_name, node2_name)] += 1
                    else:
                        E_transition[(node1_name, node2_name)] = 1
                else:
                    raise ValueError

        return E_segment, E_transition

    def get_segment_distance(self):
        """
        return the total distance of the prunned graphs' remaining segment edges
        :return:
        """
        # reverse the dict
        node_name_to_node_position = reverse_dict(self.node_name)

        total_distance = 0
        for node1, value in self.karsim_dict.items():
            for node2 in value:
                if node2[2].startswith('transition'):
                    continue
                # all segment edge are between the same chromosome nodes
                edge_type = self.segment_edge_type[(node1[0], node1[1], node2[0], node2[1])]
                # total_distance += abs(node2[1] - node1[1] + 1)
                if edge_type.startswith('telomere') or edge_type.startswith('acrocentric'):
                    continue
                else:
                    total_distance += abs(node2[1] - node1[1] + 1)

        for node1, value in self.omkar_dict.items():
            for node2 in value:
                if node2[2].startswith('transition'):
                    continue
                # all segment edge are between the same chromosome nodes
                edge_type = self.segment_edge_type[(node1[0], node1[1], node2[0], node2[1])]

                # total_distance += abs(node2[1] - node1[1] + 1)
                # TODO: recover this
                if edge_type.startswith('telomere') or edge_type.startswith('acrocentric'):
                    continue
                else:
                    total_distance += abs(node2[1] - node1[1] + 1)

        return total_distance

    def get_missed_transition_edges(self):
        """
        archived
        :return:
        """
        karsim_transition_edges = []
        for node1, value in self.karsim_dict.items():
            for node2 in value:
                if node2[2].startswith('transition'):
                    # skip reference edge
                    if node1[0] == node2[0] and node1[1] + 1 == node2[1]:
                        continue
                    else:
                        if node1 in self.node_name:
                            node1_name = self.node_name[node1]
                        else:
                            node1_name = self.source_sink_nodes[node1]
                        if node2[:2] in self.node_name:
                            node2_name = self.node_name[node2[:2]]
                        else:
                            node2_name = self.source_sink_nodes[node2[:2]]
                        karsim_transition_edges.append((node1[0], node1[1], node2[0], node2[1], f"{node1_name},{node2_name}"))

        omkar_transition_edges = []
        for node1, value in self.omkar_dict.items():
            for node2 in value:
                if node2[2].startswith('transition'):
                    # skip reference edge
                    if node1[0] == node2[0] and node1[1] + 1 == node2[1]:
                        continue
                    else:
                        if node1 in self.node_name:
                            node1_name = self.node_name[node1]
                        else:
                            node1_name = self.source_sink_nodes[node1]
                        if node2[:2] in self.node_name:
                            node2_name = self.node_name[node2[:2]]
                        else:
                            node2_name = self.source_sink_nodes[node2[:2]]
                        omkar_transition_edges.append((node1[0], node1[1], node2[0], node2[1], f"{node1_name},{node2_name}"))

        return karsim_transition_edges, omkar_transition_edges

    def get_terminally_labeled_residual_transitions(self, karsim_terminal_distance, omkar_terminal_distance):
        """
        refactor to ARCHIEVE, now uses label status and terminally labeled event
        :param terminal_distance:
        :return:
        """
        boundaries = get_prefix_suffix_forbidden_boundaries()

        def near_boundary(input_node, distance):
            input_chr = input_node[0]
            input_pos = input_node[1]
            start_boundary = boundaries[input_chr]['start']
            end_boundary = boundaries[input_chr]['end']
            if abs(start_boundary - input_pos) <= distance or abs(end_boundary - input_pos) <= distance:
                return True
            return False

        karsim_terminal_transitions = []  # dict of edge: edge label
        karsim_nonterminal_transitions = []  # same
        omkar_terminal_transitions = []   # edges; as we do not have labels for omkar side
        omkar_nonterminal_transitions = []

        for node1, node1_nbhd in self.karsim_dict.items():
            for node2 in node1_nbhd:
                if node2[2] == 'segment':
                    continue
                edge = (node1[0], node1[1], node2[0], node2[1])
                if near_boundary(node1, karsim_terminal_distance) or near_boundary(node2, karsim_terminal_distance):
                    karsim_terminal_transitions.append({edge: self.karsim_edge_label[edge]})
                else:
                    karsim_nonterminal_transitions.append({edge: self.karsim_edge_label[edge]})
        for node1, node1_nbhd in self.omkar_dict.items():
            for node2 in node1_nbhd:
                if node2[2] == 'segment':
                    continue
                edge = (node1[0], node1[1], node2[0], node2[1])
                if near_boundary(node1, omkar_terminal_distance) or near_boundary(node2, omkar_terminal_distance):
                    omkar_terminal_transitions.append(edge)
                else:
                    omkar_nonterminal_transitions.append(edge)
        return karsim_terminal_transitions, karsim_nonterminal_transitions, omkar_terminal_transitions, omkar_nonterminal_transitions

    def get_residual_edges_event_counts(self):
        type_counts = {}
        for edge, labels in self.karsim_edge_label.items():
            for label in labels:
                label_type = label.split('#')[0]
                if label_type in type_counts:
                    type_counts[label_type] += 1
                else:
                    type_counts[label_type] = 1
        return type_counts

    def get_chr_start_end_nodes(self):
        """
        find the two nodes that create a chr-breakpoint (i.e. belong to two diff. chr)
        :return: list of all chr-breakpoint nodes' name
        """
        nodes = self.node_name.keys()
        sorted_nodes = sorted(nodes, key=custom_sort_node)  # sorted by chr, then by positions -> chr always start to end sorted in the list
        terminal_nodes = [self.node_name[sorted_nodes[0]]]  # first node always the start of a chr
        node_ind = 1
        while node_ind <= len(sorted_nodes) - 3:
            current_node = sorted_nodes[node_ind]
            next_node = sorted_nodes[node_ind + 1]
            if current_node[0] != next_node[0]:
                terminal_nodes.append(self.node_name[current_node])
                terminal_nodes.append(self.node_name[next_node])
            node_ind += 2  # because of the st paired structure, only check the t against the next s
        terminal_nodes.append(self.node_name[sorted_nodes[-1]])  # last node always the end of a chr

        return terminal_nodes

    def remove_all_segment_edges(self):
        def remove_all_segment_edges_from_dict(dict_of_interest):
            node1_to_remove = []
            for node1, node2s in dict_of_interest.items():
                edges_to_remove = []
                for idx, node2 in enumerate(node2s):
                    if node2[2] == 'segment':
                        edges_to_remove.append(idx)
                dict_of_interest[node1] = [itr for idx, itr in enumerate(node2s) if idx not in edges_to_remove]
                if not dict_of_interest[node1]:
                    node1_to_remove.append(node1)
            for node1 in node1_to_remove:
                dict_of_interest.pop(node1)

        remove_all_segment_edges_from_dict(self.karsim_dict)
        remove_all_segment_edges_from_dict(self.omkar_dict)

    def remove_approximate_transition_edges(self, approximate_allowance=i_approximate_allowance):
        def mark_approximated_transition_edges(target_graph):
            if target_graph == 'karsim':
                target_dict = self.karsim_dict
            elif target_graph == 'omkar':
                target_dict = self.omkar_dict
            else:
                raise ValueError

            # mark all SMALL edges for removal
            edges_to_remove = []
            for node1, value in target_dict.items():
                for node2 in value:
                    chr1 = node1[0]
                    pos1 = node1[1]
                    chr2 = node2[0]
                    pos2 = node2[1]
                    edge_type = node2[2]

                    if edge_type == 'segment':
                        continue

                    intra_distance = intra_transition_edge_distance(chr1, pos1, chr2, pos2)
                    if intra_distance == -1:
                        # -1 encodes +inf
                        continue
                    elif is_self_edge(chr1, pos1, chr2, pos2):
                        # do not pop self-edge
                        continue
                    elif intra_distance < approximate_allowance:
                        edges_to_remove.append((chr1, pos1, chr2, pos2, edge_type))
            return edges_to_remove

        # removal
        karsim_edges_to_remove = mark_approximated_transition_edges('karsim')
        omkar_edges_to_remove = mark_approximated_transition_edges('omkar')
        removed_edge_labels = []
        for current_edge in karsim_edges_to_remove:
            distance_param = current_edge[:4]
            current_distance = intra_transition_edge_distance(*distance_param)
            self.approximated_cnv += current_distance  # refactor: archieve this summed cnv
            if current_distance > 0:
                self.karsim_n_transition_approximated += 1

            pop_edge(*current_edge, self.karsim_dict)
            if current_edge in self.karsim_edge_label:
                removed_edge_labels.append((self.karsim_edge_label[current_edge][0], current_distance))
                if len(self.karsim_edge_label[current_edge]) == 1:
                    self.karsim_edge_label.pop(*current_edge)
                else:
                    self.karsim_edge_label[current_edge].pop(0)
            else:
                removed_edge_labels.append(('ENF', current_distance))
        for current_edge in omkar_edges_to_remove:
            distance_param = current_edge[:4]
            current_distance = intra_transition_edge_distance(*distance_param)
            self.approximated_cnv += current_distance
            if current_distance > 0:
                self.omkar_n_transition_approximated += 1
            pop_edge(*current_edge, self.omkar_dict)
        return removed_edge_labels

    def node_name_to_coordinate(self, node_name):
        reverse_node_name = reverse_dict(self.node_name)
        reverse_source_sink_nodes = reverse_dict(self.source_sink_nodes)
        if node_name in reverse_node_name:
            return reverse_node_name[node_name]
        else:
            return reverse_source_sink_nodes[node_name]

    def match_transition_edges(self, matching_allowance=i_matching_allowance):
        large_value = 1000000000

        ## gather all transition edges
        _, karsim_transition_edge_dict = self.gather_edges('karsim')
        _, omkar_transition_edge_dict = self.gather_edges('omkar')
        karsim_transition_edges = []
        omkar_transition_edges = []

        for edge, multiplicity in karsim_transition_edge_dict.items():
            node1 = self.node_name_to_coordinate(edge[0])
            node2 = self.node_name_to_coordinate(edge[1])
            chr1 = node1[0]
            pos1 = node1[1]
            chr2 = node2[0]
            pos2 = node2[1]
            for itr in range(multiplicity):
                karsim_transition_edges.append((chr1, pos1, chr2, pos2))

        for edge, multiplicity in omkar_transition_edge_dict.items():
            node1 = self.node_name_to_coordinate(edge[0])
            node2 = self.node_name_to_coordinate(edge[1])
            chr1 = node1[0]
            pos1 = node1[1]
            chr2 = node2[0]
            pos2 = node2[1]
            for itr in range(multiplicity):
                omkar_transition_edges.append((chr1, pos1, chr2, pos2))

        if len(karsim_transition_edges) == 0 or len(omkar_transition_edges) == 0:
            return

        ## populate cost_matrix
        def node_sign(node_chr, node_pos):
            if (node_chr, node_pos) in self.node_name:
                return self.node_name[(node_chr, node_pos)][-1]
            else:
                if self.source_sink_nodes[(node_chr, node_pos)][0] == 'T':
                    return 's'
                else:
                    return 't'
                # return self.source_sink_nodes[(node_chr, node_pos)][0]

        cost_matrix = [[large_value for col in range(len(omkar_transition_edges))] for row in range(len(karsim_transition_edges))]
        for row_index in range(len(karsim_transition_edges)):
            for col_index in range(len(omkar_transition_edges)):
                karsim_edge = karsim_transition_edges[row_index]
                omkar_edge = omkar_transition_edges[col_index]

                # if the start/end direction do not match, transition edges are different
                karsim_node1_sign = node_sign(karsim_edge[0], karsim_edge[1])
                karsim_node2_sign = node_sign(karsim_edge[2], karsim_edge[3])
                omkar_node1_sign = node_sign(omkar_edge[0], omkar_edge[1])
                omkar_node2_sign = node_sign(omkar_edge[2], omkar_edge[3])
                if karsim_node1_sign != omkar_node1_sign or karsim_node2_sign != omkar_node2_sign:
                    continue

                distance = inter_transition_edge_distance(*karsim_edge, *omkar_edge)
                if distance == -1:
                    # -1 signals no match
                    continue

                if distance < 2 * matching_allowance:
                    # only update cost if within threshold
                    cost_matrix[row_index][col_index] = distance

        ## bipartite matching
        np_cost_matrix = np.array(cost_matrix)
        karsim_assignment, omkar_assignment = linear_sum_assignment(np_cost_matrix)

        ## remove all approximate matching
        approximate_matching_labels = []  # (label, distance)
        for ind in range(len(karsim_assignment)):
            karsim_ind = karsim_assignment[ind]
            omkar_ind = omkar_assignment[ind]
            current_distance = cost_matrix[karsim_ind][omkar_ind]

            if current_distance < large_value:
                ## an approximate matching is found, pop both edges in the matching
                self.approximated_cnv += current_distance
                karsim_edge = karsim_transition_edges[karsim_ind]
                omkar_edge = omkar_transition_edges[omkar_ind]
                ## record distances in unit of each coordinate, for self-edge, record only one distance
                # for a pair of edge to be matched, it is guarenteed to be in the same orientation and same Chr-pairings
                if karsim_edge[0:2] == karsim_edge[2:] and omkar_edge[0:2] == omkar_edge[2:]:
                    # self-edge
                    distances_by_each_coordinate = [abs(karsim_edge[1] - omkar_edge[1])]
                else:
                    distances_by_each_coordinate = [abs(karsim_edge[1] - omkar_edge[1]), abs(karsim_edge[3] - omkar_edge[3])]

                ## document edge removed, all edges now are significant and should have a label
                if karsim_edge in self.karsim_edge_label:
                    approximate_matching_labels.append({'event_id': self.karsim_edge_label[karsim_edge][0],
                                                        'distances': distances_by_each_coordinate,
                                                        'karsim_edge': karsim_edge,
                                                        'omkar_edge': omkar_edge})
                    if len(self.karsim_edge_label[karsim_edge]) == 1:
                        self.karsim_edge_label.pop(karsim_edge)
                    else:
                        # (very rare) this edge has multiple event-id label, we can't tell which label was matched, so we just pop one of the label
                        self.karsim_edge_label[karsim_edge].pop(0)
                else:
                    approximate_matching_labels.append({'event_id': "ENF",
                                                        'distances': distances_by_each_coordinate,
                                                        'karsim_edge': karsim_edge,
                                                        'omkar_edge': omkar_edge})

                ## karsim edge removal
                c_dict_entry = self.karsim_dict[(karsim_edge[0], karsim_edge[1])]
                if (karsim_edge[2], karsim_edge[3], 'transition') in c_dict_entry:
                    c_dict_entry.remove((karsim_edge[2], karsim_edge[3], 'transition'))
                elif (karsim_edge[2], karsim_edge[3], 'transition_source') in c_dict_entry:
                    c_dict_entry.remove((karsim_edge[2], karsim_edge[3], 'transition_source'))
                elif (karsim_edge[2], karsim_edge[3], 'transition_sink') in c_dict_entry:
                    c_dict_entry.remove((karsim_edge[2], karsim_edge[3], 'transition_sink'))
                else:
                    raise RuntimeError()
                ## omkar edge removal
                c_dict_entry = self.omkar_dict[(omkar_edge[0], omkar_edge[1])]
                if (omkar_edge[2], omkar_edge[3], 'transition') in c_dict_entry:
                    c_dict_entry.remove((omkar_edge[2], omkar_edge[3], 'transition'))
                elif (omkar_edge[2], omkar_edge[3], 'transition_source') in c_dict_entry:
                    c_dict_entry.remove((omkar_edge[2], omkar_edge[3], 'transition_source'))
                elif (omkar_edge[2], omkar_edge[3], 'transition_sink') in c_dict_entry:
                    c_dict_entry.remove((omkar_edge[2], omkar_edge[3], 'transition_sink'))
                else:
                    raise RuntimeError()

                self.karsim_n_transition_approximated += 1
                self.omkar_n_transition_approximated += 1
        return approximate_matching_labels

    def node_is_small(self, internal_node_name, size_threshold):
        """
        only works with internal nodes (not source/sink)
        require running generate_node_name_reverse_dict prior to calling this
        :param size_threshold: to be labeled as small
        :param internal_node_name:
        :return:
        """
        node_base_name = internal_node_name[:-1]
        node_side = internal_node_name[-1]
        if node_side == 't':
            paired_node_name = node_base_name + 's'
        elif node_side == 's':
            paired_node_name = node_base_name + 't'
        else:
            raise RuntimeError("illegal node suffix")
        node1_pos = self.node_name_reverse_dict[internal_node_name][1]
        node2_pos = self.node_name_reverse_dict[paired_node_name][1]
        if abs(node1_pos - node2_pos) <= size_threshold:
            return True
        else:
            return False

    def remove_forbidden_nodes(self, forbidden_region_file):
        """
        remove all nodes associated with forbidden regions, and then remove all edges associated with these nodes
        :param forbidden_region_file:
        :return:
        """
        def remove_node_and_associated_edge(input_node):
            self.node_name.pop(input_node)
            # remove out going edges
            if input_node in self.karsim_dict:
                self.karsim_dict.pop(input_node)
            if input_node in self.omkar_dict:
                self.omkar_dict.pop(input_node)
            # remove incoming edges
            for node1, node2_list in self.karsim_dict.items():
                node2_to_remove = []
                for idx, node2 in enumerate(node2_list):
                    if node2[0] == input_node[0] and node2[1] == input_node[1]:
                        node2_to_remove.append(idx)
                new_node2_list = [value for i, value in enumerate(node2_list) if i not in node2_to_remove]
                self.karsim_dict[node1] = new_node2_list
            for node1, node2_list in self.omkar_dict.items():
                node2_to_remove = []
                for idx, node2 in enumerate(node2_list):
                    if node2[0] == input_node[0] and node2[1] == input_node[1]:
                        node2_to_remove.append(idx)
                new_node2_list = [value for i, value in enumerate(node2_list) if i not in node2_to_remove]
                self.omkar_dict[node1] = new_node2_list

        forbidden_segments = read_forbidden_regions(forbidden_region_file).segments
        for segment in forbidden_segments:
            start_node = (segment.chr_name, segment.start)
            end_node = (segment.chr_name, segment.end)
            if start_node in self.node_name:
                remove_node_and_associated_edge(start_node)
            if end_node in self.node_name:
                remove_node_and_associated_edge(end_node)

    def visualize_graph(self, output_prefix, merged=False):
        from matplotlib import rcParams
        rcParams['pdf.fonttype'] = 42

        # create sorted nodes (Endpoints)

        nodes = list(self.node_name.keys()) + list(self.source_sink_nodes.keys())
        sorted_nodes = sorted(nodes, key=custom_sort_node)
        V = []
        for node in sorted_nodes:
            if node in self.node_name:
                V.append(self.node_name[node])
            elif node in self.source_sink_nodes:
                V.append(self.source_sink_nodes[node])
            else:
                raise RuntimeError('node not in the total node dict')

        E_karsim_segment, E_karsim_transition = self.gather_edges('karsim')
        E_omkar_segment, E_omkar_transition = self.gather_edges('omkar')

        def translate_segment_edge_type_to_edge_name():
            new_dict = {}
            for edge, edge_type in self.segment_edge_type.items():
                if (edge[0], edge[1]) in self.node_name and (edge[2], edge[3]) in self.node_name:
                    node1_name = self.node_name[(edge[0], edge[1])]
                    node2_name = self.node_name[(edge[2], edge[3])]
                    new_dict[(node1_name, node2_name)] = edge_type
            return new_dict

        def iterative_add_edge(edges_dict: {(str, str): int}, group_color, graph, forbidden_segment_edge_labels=False):
            """
            if forbidden_segment_edge_labels=True, prepare "translated_segment_edge_type"
            :param edges_dict:
            :param group_color:
            :param graph:
            :param forbidden_segment_edge_labels:
            :return:
            """
            for edge, edge_weight in edges_dict.items():
                edge_color = group_color
                if forbidden_segment_edge_labels:
                    this_edge_type = translated_segment_edge_type[(edge[0]), (edge[1])]
                    if this_edge_type.startswith('telomere') or this_edge_type.startswith('acrocentric'):
                        edge_color = 'blue'

                if graph.has_edge(edge[0], edge[1]):
                    # FIXME: right now, if different types of the same edge exist, the multiplicity of the latters will be added to the first added type
                    graph[edge[0]][edge[1]]['weight'] += edge_weight
                else:
                    graph.add_edge(edge[0], edge[1], color=edge_color, weight=edge_weight)

        def iterative_add_edge_trace(edges_dict: {(str, str): int},
                                     trace_dict: {(str, str): (float, float)},
                                     peak_height_nonadjacent,
                                     peak_height_adjacent,
                                     adjacent_x_dist):
            """
            :param edges_dict:
            :param trace_dict:
            :param peak_height_nonadjacent: for nodes not right next to each other
            :param peak_height_adjacent: for nodes right next to each other
            :param adjacent_x_dist: nodes x-dist if they are adjacent
            :return:
            """
            for edge in edges_dict:
                node1_x = V_pos[edge[0]][0]
                node2_x = V_pos[edge[1]][0]

                if edge in trace_dict:
                    ## same edge of multiple type will be colored by the first type's color
                    continue

                if edge[0] == edge[1]:
                    # self edge
                    trace_dict[edge] = generate_circle(node1_x, peak_height_nonadjacent, 1)
                else:
                    # regular edge
                    if abs(node2_x - node1_x) <= adjacent_x_dist + 0.001:
                        peak_x = round((node1_x + node2_x) / 2, 6)
                        trace_dict[edge] = generate_parabola(node1_x, node2_x, peak_x, peak_height_adjacent)
                    else:
                        peak_x = round((node1_x + node2_x) / 2, 6)
                        trace_dict[edge] = generate_parabola(node1_x, node2_x, peak_x, peak_height_nonadjacent)

        G_karsim = nx.DiGraph()
        G_karsim.add_nodes_from(V)

        G_omkar = nx.DiGraph()
        G_omkar.add_nodes_from(V)

        # add all edges
        translated_segment_edge_type = translate_segment_edge_type_to_edge_name()
        iterative_add_edge(E_karsim_segment, 'black', G_karsim, forbidden_segment_edge_labels=True)
        iterative_add_edge(E_karsim_transition, 'red', G_karsim)
        iterative_add_edge(E_omkar_segment, 'black', G_omkar, forbidden_segment_edge_labels=True)
        iterative_add_edge(E_omkar_transition, 'red', G_omkar)

        def generate_circular_coordinates(n, radius=0.5, center=(0.5, 0.5)):
            cx, cy = center
            coordinates = []

            for k in range(n):
                angle = (k / n) * 2 * math.pi
                x = cx + radius * math.cos(angle)
                y = cy + radius * math.sin(angle)
                coordinates.append((x, y))

            return coordinates

        def generate_linear_coordinates(n, left_boundary=0, right_boundary=5, fixed_y=0.5):
            coordinates = []
            for i in range(n):
                x = 0 + (i / (n - 1)) * (right_boundary - left_boundary)
                y = fixed_y
                coordinates.append((x, y))

            return coordinates

        uniform_dist = 1 / 6
        # graph_width = uniform_dist * 200
        # width_multiplier = 2.5
        # height_multiplier = 1

        def generate_uniform_linear_coordinates(n, fixed_y=0.5, fixed_distance=uniform_dist):
            coordinates = []
            for i in range(n):
                x = round(0 + i * fixed_distance, 4)
                y = fixed_y
                coordinates.append((x, y))

            return coordinates

        ## node parameters: same for the two graphs (nodes are always the same)
        V_pos = {}
        # cor = generate_circular_coordinates(len(V))
        cor = generate_uniform_linear_coordinates(len(V))
        for node_itr_ind in range(len(V)):
            V_pos[V[node_itr_ind]] = cor[node_itr_ind]

        V_colors = {}
        internal_nodes = list(self.node_name.values())
        source_sink_nodes = list(self.source_sink_nodes.values())
        # terminal_nodes = self.get_chr_start_end_nodes()
        self.generate_node_name_reverse_dict()
        for node_itr in internal_nodes:
            if self.node_is_small(node_itr, size_threshold=i_approximate_allowance):
                V_colors[node_itr] = 'cyan'
            else:
                V_colors[node_itr] = 'white'
        for node_itr in source_sink_nodes:
            V_colors[node_itr] = 'orange'

        ## edges parameters: diff. for the two graphs
        karsim_E_weights = {}
        karsim_E_colors = nx.get_edge_attributes(G_karsim, 'color')
        for edge_itr in E_karsim_segment:
            karsim_E_weights[edge_itr] = E_karsim_segment[edge_itr]
        for edge_itr in E_karsim_transition:
            karsim_E_weights[edge_itr] = E_karsim_transition[edge_itr]

        omkar_E_weights = {}
        omkar_E_colors = nx.get_edge_attributes(G_omkar, 'color')
        for edge_itr in E_omkar_segment:
            omkar_E_weights[edge_itr] = E_omkar_segment[edge_itr]
        for edge_itr in E_omkar_transition:
            omkar_E_weights[edge_itr] = E_omkar_transition[edge_itr]

        ## plotting
        if not merged:
            fixed_height = 6
            karsim_E_pos = {}
            iterative_add_edge_trace(E_karsim_segment, karsim_E_pos, 1, 0.58, uniform_dist)
            iterative_add_edge_trace(E_karsim_transition, karsim_E_pos, 1, 0.58, uniform_dist)

            omkar_E_pos = {}
            iterative_add_edge_trace(E_omkar_segment, omkar_E_pos, 0, 0.42, uniform_dist)
            iterative_add_edge_trace(E_omkar_transition, omkar_E_pos, 0, 0.42, uniform_dist)

            graph_width = len(V_pos)
            plt.figure(figsize=(graph_width, fixed_height))
            plot_karsim = ng.InteractiveGraph(G_karsim,
                                              node_color=V_colors,
                                              node_layout=V_pos,
                                              node_labels=True,
                                              edge_color=karsim_E_colors,
                                              edge_layout=karsim_E_pos,
                                              edge_labels=karsim_E_weights,
                                              arrows=True,
                                              node_size=6,
                                              node_label_offset=0.001,
                                              node_label_font_dict=dict(size=20),
                                              edge_label_fontdict=dict(size=15),
                                              scale=(graph_width, 1))
            plt.savefig(output_prefix + '.karsim.graph.pdf', dpi=500)

            plt.figure(figsize=(graph_width, fixed_height))
            plot_omkar = ng.InteractiveGraph(G_omkar,
                                             node_color=V_colors,
                                             node_layout=V_pos,
                                             node_labels=True,
                                             edge_color=omkar_E_colors,
                                             edge_layout=omkar_E_pos,
                                             edge_labels=omkar_E_weights,
                                             arrows=True,
                                             node_size=6,
                                             node_label_offset=0.001,
                                             node_label_font_dict=dict(size=20),
                                             edge_label_fontdict=dict(size=15),
                                             scale=(graph_width, 1))
            plt.savefig(output_prefix + '.omkar.graph.pdf', dpi=500)
        else:
            G_merged = nx.DiGraph()
            G_merged.add_nodes_from(V)

            iterative_add_edge(E_karsim_segment, 'black', G_merged, forbidden_segment_edge_labels=True)
            iterative_add_edge(E_karsim_transition, 'red', G_merged)
            iterative_add_edge(E_omkar_segment, 'black', G_merged, forbidden_segment_edge_labels=True)
            iterative_add_edge(E_omkar_transition, 'red', G_merged)

            merged_E_pos = {}
            adjacent_trace_height = 0.08
            adjacent_transition_multiplier = 2
            # adjacent_transition_multiplier = 1
            nonadjacent_trace_height = min(max(50, len(V_pos)), 150) * 0.01
            # nonadjacent_trace_height = 0.25
            iterative_add_edge_trace(E_karsim_segment, merged_E_pos, 0.5 + nonadjacent_trace_height, 0.5 + adjacent_trace_height, uniform_dist)
            iterative_add_edge_trace(E_karsim_transition, merged_E_pos, 0.5 + nonadjacent_trace_height, 0.5 + adjacent_trace_height * adjacent_transition_multiplier, uniform_dist)
            iterative_add_edge_trace(E_omkar_segment, merged_E_pos, 0.5 - nonadjacent_trace_height, 0.5 - adjacent_trace_height, uniform_dist)
            iterative_add_edge_trace(E_omkar_transition, merged_E_pos, 0.5 - nonadjacent_trace_height, 0.5 - adjacent_trace_height * adjacent_transition_multiplier, uniform_dist)

            merged_E_colors = {**karsim_E_colors, **omkar_E_colors}
            merged_E_weights = {**karsim_E_weights, **omkar_E_weights}



            graph_width = len(V_pos)
            fixed_height = min(max(60, len(V_pos)), 150) * 0.1  # to scale height according to the trace height
            plt.figure(figsize=(graph_width, fixed_height * 1.8))  # merged graph is going to be twice as high
            plot_merged = ng.InteractiveGraph(G_merged,
                                              node_color=V_colors,
                                              node_layout=V_pos,
                                              node_labels=True,
                                              edge_color=merged_E_colors,
                                              edge_layout=merged_E_pos,
                                              edge_labels=merged_E_weights,
                                              arrows=True,
                                              node_size=6,
                                              node_label_fontdict={'fontsize': 'large'},
                                              node_label_offset=0.001,
                                              node_label_font_dict=dict(size=10),
                                              edge_label_fontdict=dict(size=9),
                                              scale=(graph_width, 1))
            plt.savefig(output_prefix + '.merged.graph.pdf', dpi=500)

    def translate_transition_edges(self):
        """
        for debug use
        :return:
        """
        for node1, node1_nbhd in self.karsim_dict.items():
            new_node1 = (node1[0], node1[1], self.node_name[node1])  # chr, pos, node_name
            new_node1_nbhd = []
            for node2 in node1_nbhd:
                if not node2[2].startswith('transition'):
                    continue
                if node2[2] == 'transition':
                    new_node1_nbhd.append((node2[0], node2[1], self.node_name[node2[:2]]))
                else:
                    new_node1_nbhd.append((node2[0], node2[1], self.source_sink_nodes[node2[:2]]))
            if new_node1_nbhd:
                self.translated_karsim_dict[new_node1] = new_node1_nbhd
        for node1, node1_nbhd in self.omkar_dict.items():
            new_node1 = (node1[0], node1[1], self.node_name[node1])
            new_node1_nbhd = []
            for node2 in node1_nbhd:
                if not node2[2].startswith('transition'):
                    continue
                if node2[2] == 'transition':
                    new_node1_nbhd.append((node2[0], node2[1], self.node_name[node2[:2]]))
                else:
                    new_node1_nbhd.append((node2[0], node2[1], self.source_sink_nodes[node2[:2]]))
            if new_node1_nbhd:
                self.translated_omkar_dict[new_node1] = new_node1_nbhd

    def report_chr_in_graph(self):
        chroms = set()
        for node1, node1_nbhd in self.karsim_dict.items():
            chroms.add(node1[0])
            for node2 in node1_nbhd:
                chroms.add(node2[0])
        for node1, node1_nbhd in self.omkar_dict.items():
            chroms.add(node1[0])
            for node2 in node1_nbhd:
                chroms.add(node2[0])
        return chroms

def tally_graph_edge_labels(input_graph):
    """
    report SV by SV-id, their multiplicity, and the current distance for the SV
    :param input_graph:
    :return:
    """
    output_dict = {}
    for edge, labels in input_graph.karsim_edge_label.items():
        for label in labels:
            if label not in output_dict:
                output_dict[label] = {'count': 1, 'distances': []}  # multiplicity, matched_distances_by_edge ([int], each entry is the coordinate distance)
            else:
                output_dict[label]['count'] += 1
    return output_dict

def label_status_log_initial_count(input_labels):
    for event_id, status_entry in input_labels.items():
        status_entry['initial count'] = status_entry['count']

def label_status_populate_edge_fields(input_labels):
    for event_id in input_labels:
        input_labels[event_id]['karsim_edge'] = []
        input_labels[event_id]['omkar_edge'] = []

def update_label_status(input_labels, removed_labels):
    if removed_labels is None:
        return
    for entry in removed_labels:
        if entry['event_id'] != 'ENF':
            input_labels[entry['event_id']]['count'] -= 1
            input_labels[entry['event_id']]['distances'].append(entry['distances'])
            input_labels[entry['event_id']]['karsim_edge'].append(entry['karsim_edge'])
            input_labels[entry['event_id']]['omkar_edge'].append(entry['omkar_edge'])
        else:
            if 'ENF' not in input_labels:
                input_labels['ENF'] = {'count': -1, 'distances': []}
            else:
                input_labels['ENF']['count'] += -1
                input_labels['ENF']['distances'].append(entry['distances'])
                input_labels['ENF'].append(entry['karsim_edge'])
                input_labels['ENF'].append(entry['omkar_edge'])

def update_terminal_label_status(input_labels, removed_labels):
    ## this version allows the removed label to not be a subgroup of the input_labels
    if removed_labels is None:
        return
    for entry in removed_labels:
        c_label = entry[0]
        c_distance = entry[1]
        if c_label in input_labels:
            input_labels[c_label]['count'] -= 1
            input_labels[c_label]['distances'].append(c_distance)

def deep_copy_status(input_label_status):
    output_label_status = {}
    for event_id, info in input_label_status.items():
        output_label_status[event_id] = {'count': info['count'], 'distances': copy.deepcopy(info['distances'])}
    return output_label_status

def split_event_status_into_terminal_nonterminal(terminal_event_ids, input_event_status):
    T_event_status = {}
    nT_event_status = {}
    for event_id, status in input_event_status.items():
        if event_id in terminal_event_ids:
            T_event_status[event_id] = status
        else:
            nT_event_status[event_id] = status
    return T_event_status, nT_event_status

def split_event_status_into_caught_uncaught_events(input_event_status):
    caught_event_status = {}
    uncaught_event_status = {}
    for event_id, status in input_event_status.items():
        if status['count'] == 0:
            caught_event_status[event_id] = status
        else:
            uncaught_event_status[event_id] = status
    return caught_event_status, uncaught_event_status

def event_status_get_caught_uncaught_counts(input_event_status):
    caught_events = 0
    uncaught_events = 0
    caught_edges = 0
    uncaught_edges = 0
    for event_id, status in input_event_status.items():
        if status['count'] == 0:
            caught_events += 1
        else:
            uncaught_events += 1
        caught_edges += status['initial count'] - status['count']
        uncaught_edges += status['count']
    return caught_events, uncaught_events, caught_edges, uncaught_edges

def form_graph_from_cluster(cluster_file):
    """
    for graph without the prefix/suffix forbidden regions,
    and then add source/sink nodes for prefix/suffix transitions into the middle nonforbidden regions
    :param cluster_file:
    :return:
    """
    index_to_segment, karsim_path_list, omkar_path_list, labeled_edges = read_cluster_file(cluster_file)
    graph = Graph()
    graph.edges_of_interest = labeled_edges

    # filter out prefix/suffix forbidden segments; forbidden segment's breakpoints are already in-place
    for path in karsim_path_list:
        path.drop_forbidden_region_segments()
    for path in omkar_path_list:
        path.drop_forbidden_region_segments()

    def add_segment_edge(path_list, target_graph):
        for path_itr in path_list:
            for segment in path_itr.linear_path.segments:
                graph.add_segment_edge(segment, target_graph)

    def add_transition_edge(path_list, target_graph):
        for path_itr in path_list:
            for segment_ind in range(len(path_itr.linear_path.segments) - 1):
                current_segment = path_itr.linear_path.segments[segment_ind]
                next_segment = path_itr.linear_path.segments[segment_ind + 1]
                graph.add_transition_edge(current_segment, next_segment, target_graph)

    # segment edge
    add_segment_edge(karsim_path_list, 'karsim')
    add_segment_edge(omkar_path_list, 'omkar')
    # transition edge
    add_transition_edge(karsim_path_list, 'karsim')
    add_transition_edge(omkar_path_list, 'omkar')

    # add source and sink nodes, and their transitions
    boundaries = get_prefix_suffix_forbidden_boundaries()
    source_nodes = set()
    sink_nodes = set()

    def add_source_and_sink_transition_edge(path_list, target_graph):
        for path_itr in path_list:
            start_seg = path_itr.linear_path.segments[0]
            end_seg = path_itr.linear_path.segments[-1]

            source_chr = start_seg.chr_name
            source_pos = boundaries[source_chr]['start']
            sink_chr = end_seg.chr_name
            sink_pos = boundaries[sink_chr]['end']
            source_nodes.add((source_chr, source_pos))
            sink_nodes.add((sink_chr, sink_pos))

            graph.add_source_transition_edge(source_pos, start_seg, target_graph)
            graph.add_sink_transition_edge(sink_pos, end_seg, target_graph)

    add_source_and_sink_transition_edge(karsim_path_list, 'karsim')
    add_source_and_sink_transition_edge(omkar_path_list, 'omkar')

    # add source and sink nodes
    source_nodes = list(source_nodes)
    sink_nodes = list(sink_nodes)
    source_nodes = sorted(source_nodes, key=custom_sort_node)
    sink_nodes = sorted(sink_nodes, key=custom_sort_node)
    for idx, node in enumerate(source_nodes):
        graph.source_sink_nodes[node] = f"S{idx}"
    for idx, node in enumerate(sink_nodes):
        graph.source_sink_nodes[node] = f"T{idx}"

    return graph

def chry_in_genome(input_omkar_path_list):
    """
    heuristics used for assigning correct WT sex
    :param input_omkar_path_list:
    :return: bool
    """
    for path in input_omkar_path_list:
        if "ChrY" in path.path_chr:
            return True
    return False

def compare_aneuplodies(cluster_file):
    ## for cluster where the karsim path_list has a devidation from diploid, we report TP/FN from the OMKar's output
    ## for cluster where the karsim path_list has no deviation, we report TN/FP from the OMKar's output
    ## currently ignore aneuploidies on X (for simplicity of implementation)
    # i.e. for 1-22, WT=2; X, WT={1,2}; Y, WT={0,1}
    index_to_segment, karsim_path_list, omkar_path_list, labeled_edges = read_cluster_file(cluster_file)

    # bin paths to get the Chr-origins
    rotate_and_bin_path(karsim_path_list)
    rotate_and_bin_path(omkar_path_list)

    karsim_chr_origin = {}
    for path in karsim_path_list:
        chr_origin = path.path_chr.split(": ")[-1].replace('Chr', '')
        if chr_origin in karsim_chr_origin:
            karsim_chr_origin[chr_origin] += 1
        else:
            karsim_chr_origin[chr_origin] = 1

    omkar_chr_origin = {}
    for path in omkar_path_list:
        chr_origin = path.path_chr.split(": ")[-1].replace('Chr', '')
        if chr_origin in omkar_chr_origin:
            omkar_chr_origin[chr_origin] += 1
        else:
            omkar_chr_origin[chr_origin] = 1

    wt_counts = {str(i): [2] for i in range(1, 23)}
    wt_counts['X'] = [1, 2]
    wt_counts['Y'] = [0, 1]
    output = {'TN': 0, 'TP': 0, 'FP': 0, 'FN': 0}
    is_aneuploidy_positive = False
    for chr_origin, karsim_count in karsim_chr_origin.items():
        acceptable_wt_count = wt_counts[chr_origin]
        if karsim_count not in acceptable_wt_count:
            is_aneuploidy_positive = True
            break

    omkar_count = sum(omkar_chr_origin.values())
    karsim_count = sum(karsim_chr_origin.values())
    if is_aneuploidy_positive:
        # TP/FN
        if omkar_count != karsim_count:
            output['FN'] += abs(karsim_count - omkar_count)
        else:
            output['TP'] += 1
    else:
        # TN/FP
        if omkar_count != karsim_count:
            output['FP'] += abs(karsim_count - omkar_count)
        else:
            output['TN'] += 1

    return output


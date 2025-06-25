from Karsimulator_Start_Genome import *
from KarInterpreter import *
from KarUtils import *
import os


def read_KarSimulator_output_to_path(KarSimulator_output_file, forbidden_region_file):
    index_dict, path_list, event_histories = read_KarSimulator_output(KarSimulator_output_file, forbidden_region_file)
    label_path_with_forbidden_regions(path_list, forbidden_region_file)
    return index_dict, path_list


def generate_history_SV_edge_labels(event_histories, sv_edge_list, output_file_path):
    with open(output_file_path, 'w') as fp_write:
        for idx, hist in enumerate(event_histories):
            fp_write.write(str(hist) + ": " + str(sv_edge_list[idx]) + '\n')


def read_KarSimulator_output(KT_file, masking_file):
    genome = generate_genome_from_KT(KT_file)
    index_dict = genome.segment_indexing()
    t2_segments = get_t2_segments(masking_file)  # required as t2_segments start is unknown from KT file alone
    path_list = []

    for chromosome_itr in genome:
        if chromosome_itr.deleted:
            continue

        ## T1
        t1_segment = Segment(chromosome_itr.name[:-1], 0, chromosome_itr.t1_len - 1, "telomere1")
        segment_list = [t1_segment]

        ## P-arm
        for segment_itr in chromosome_itr.p_arm.segments:
            segment_itr.segment_type = "arm_region"
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        ## CEN
        cen_segment = chromosome_itr.centromere.segments[0].duplicate()
        cen_segment.segment_type = "centromere"
        cen_segment.kt_index = get_kt_index(index_dict, cen_segment)
        segment_list.append(cen_segment)

        ## Q-arm
        for segment_itr in chromosome_itr.q_arm.segments:
            segment_itr.segment_type = "arm_region"
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        ## T2
        current_t2_segment = None
        for find_t2_segment_itr in t2_segments.segments:
            if find_t2_segment_itr.chr_name == chromosome_itr.name[:-1]:
                current_t2_segment = find_t2_segment_itr
                break
        segment_list.append(current_t2_segment)
        path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))

    return index_dict, path_list, genome.histories


def typed_segments_to_indexed_segments(typed_segment_list, index_dict):
    indexed_segment_list = []
    for seg in typed_segment_list:
        temp_seg = seg.duplicate()
        if not temp_seg.direction:
            temp_seg.invert()
            direction = '-'
        else:
            direction = '+'
        c_seg_index = str(index_dict[temp_seg]) + direction
        indexed_segment_list.append(c_seg_index)
    return indexed_segment_list


def indexed_segments_to_typed_segments(indexed_segment_list, index_to_segment_dict):
    typed_segment_list = []
    for seg in indexed_segment_list:
        direction = seg[-1]
        temp_seg = index_to_segment_dict[seg[:-1]].duplicate()
        if direction == '-':
            temp_seg.invert()
        typed_segment_list.append(temp_seg)
    return typed_segment_list


def form_mt_wt_indexed_paths_for_interpreter(path, index_dict):
    size_dict = {}
    chr_bin = path.path_chr
    chr_segments = []
    for seg, seg_index in index_dict.items():
        if seg.chr_name == chr_bin:
            chr_segments.append(seg)
    chr_segments = sorted(chr_segments)
    chr_segments_indexed = []
    for seg in chr_segments:
        chr_segments_indexed.append(index_dict[seg])
        size_dict[str(index_dict[seg])] = len(seg)

    wt_segments = []
    mt_segments = []
    for seg_index in chr_segments_indexed:
        wt_segments.append(seg_index + '+')

    for seg in path.linear_path.segments:
        print(seg)
        temp_seg = seg.duplicate()
        if not temp_seg.direction():
            temp_seg.invert()
            direction = '-'
        else:
            direction = '+'
        c_seg_index = str(index_dict[temp_seg]) + direction
        mt_segments.append(c_seg_index)
        size_dict[str(index_dict[temp_seg])] = len(temp_seg)

    return wt_segments, mt_segments, size_dict


def conglomerate_segments_to_dict(previous_segment_to_index_dict, path_list):
    all_segments = set(previous_segment_to_index_dict.keys())
    for path in path_list:
        for segment in path.linear_path.segments:
            temp_segment = segment.duplicate()
            if not temp_segment.direction():
                temp_segment.invert()
            all_segments.add(temp_segment)
    all_segments = sorted(list(all_segments))
    c_idx = 1
    new_dict = {}
    for seg in all_segments:
        new_dict[seg] = str(c_idx)
        c_idx += 1
    return new_dict


def label_event_sv_edge(old_segment_to_index_dict, path_list, history_list):
    def find_path(path_name):
        for path_idx, path in enumerate(path_list):
            if path.path_name == path_name:
                return path_idx
        return None  # likely chromosomal deletion

    def find_seg_sublist(superlist_segments, sublist_segments):
        superlist_hash = [hash(seg) for seg in superlist_segments]
        sublist_hash = [hash(seg) for seg in sublist_segments]
        idx_found = -1
        for i in range(len(superlist_hash) - len(sublist_hash) + 1):
            if superlist_hash[i: i+len(sublist_hash)] == sublist_hash:
                if idx_found != -1:
                    # duplicate entry exists
                    return None
                else:
                    idx_found = i
        if idx_found == -1:
            return None
        else:
            return idx_found

    def invert_typed_segment_list(input_segs):
        inverted_segs = []
        for seg in reversed(input_segs):
            new_seg = seg.duplicate()
            new_seg.invert()
            inverted_segs.append(new_seg)
        return inverted_segs

    def invert_indexed_segment_list(input_segs):
        inverted_segs = []
        for seg in reversed(input_segs):
            if seg[-1] == '-':
                new_seg = seg[:-1] + '+'
            else:
                new_seg = seg[:-1] + '-'
            inverted_segs.append(new_seg)
        return inverted_segs

    old_index_to_segment_dict = reverse_dict(old_segment_to_index_dict)
    new_segment_to_index_dict = conglomerate_segments_to_dict(old_segment_to_index_dict, path_list)
    new_index_to_segment_dict = reverse_dict(new_segment_to_index_dict)

    sv_edges_list = []
    for hist_entry in history_list:
        event_type = hist_entry[0]
        c_sv_edges = []
        if event_type.startswith('chromosomal'):
            sv_edges_list.append(['chromosomal events'])
            continue
        c_path_from_idx = find_path(hist_entry[1])
        c_path_to_idx = find_path(hist_entry[2])
        if event_type == 'inversion':
            if c_path_from_idx is None:
                # Chr-deletion
                sv_edges_list.append(['chromosome deleted'])
                continue
            else:
                c_path_from = path_list[c_path_from_idx]
            event_segments = hist_entry[3]
            event_typed_segments = indexed_segments_to_typed_segments(event_segments, old_index_to_segment_dict)
            inverted_event_segments = invert_typed_segment_list(event_typed_segments)
            sublist_idx = find_seg_sublist(c_path_from.linear_path.segments, inverted_event_segments)
            if sublist_idx is not None:
                seg1 = c_path_from.linear_path.segments[sublist_idx - 1]
                seg2 = c_path_from.linear_path.segments[sublist_idx]
                seg3 = c_path_from.linear_path.segments[sublist_idx + len(inverted_event_segments) - 1]
                seg4 = c_path_from.linear_path.segments[sublist_idx + len(inverted_event_segments)]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'tandem duplication':
            if c_path_from_idx is None:
                sv_edges_list.append(['chromosome deleted'])
                continue
            else:
                c_path_from = path_list[c_path_from_idx]
            indexed_event_segments = hist_entry[3] + hist_entry[3]
            typed_event_segments = indexed_segments_to_typed_segments(indexed_event_segments, old_index_to_segment_dict)
            sublist_idx = find_seg_sublist(c_path_from.linear_path.segments, typed_event_segments)
            if sublist_idx is not None:
                seg1 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) // 2 - 1]
                seg2 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) // 2]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
        elif event_type == 'left duplication inversion':
            if c_path_from_idx is None:
                sv_edges_list.append(['chromosome deleted'])
                continue
            else:
                c_path_from = path_list[c_path_from_idx]
            indexed_event_segments = invert_indexed_segment_list(hist_entry[3]) + hist_entry[3]
            typed_event_segments = indexed_segments_to_typed_segments(indexed_event_segments, old_index_to_segment_dict)
            sublist_idx = find_seg_sublist(c_path_from.linear_path.segments, typed_event_segments)
            if sublist_idx is not None:
                seg1 = c_path_from.linear_path.segments[sublist_idx - 1]
                seg2 = c_path_from.linear_path.segments[sublist_idx]
                seg3 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) // 2 - 1]
                seg4 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) // 2]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'right duplication inversion':
            if c_path_from_idx is None:
                sv_edges_list.append(['chromosome deleted'])
                continue
            else:
                c_path_from = path_list[c_path_from_idx]
            indexed_event_segments = hist_entry[3] + invert_indexed_segment_list(hist_entry[3])
            typed_event_segments = indexed_segments_to_typed_segments(indexed_event_segments, old_index_to_segment_dict)
            sublist_idx = find_seg_sublist(c_path_from.linear_path.segments, typed_event_segments)
            if sublist_idx is not None:
                seg1 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) // 2 - 1]
                seg2 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) // 2]
                seg3 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments) - 1]
                seg4 = c_path_from.linear_path.segments[sublist_idx + len(indexed_event_segments)]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'balanced reciprocal translocation':
            # requirement: balanced reciprocal translocation always appear in pairs
            if c_path_to_idx is None:
                sv_edges_list.append(['chromosome deleted'])
                continue
            else:
                c_path_to = path_list[c_path_to_idx]
            indexed_event_segments = hist_entry[3]
            typed_event_segments = indexed_segments_to_typed_segments(indexed_event_segments, old_index_to_segment_dict)
            sublist_idx = find_seg_sublist(c_path_to.linear_path.segments, typed_event_segments)
            if sublist_idx is not None:
                seg1 = c_path_to.linear_path.segments[sublist_idx - 1]
                seg2 = c_path_to.linear_path.segments[sublist_idx]
                seg3 = c_path_to.linear_path.segments[sublist_idx + len(indexed_event_segments) - 1]
                seg4 = c_path_to.linear_path.segments[sublist_idx + len(indexed_event_segments)]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'deletion':
            if c_path_from_idx is None:
                sv_edges_list.append(['chromosome deleted'])
                continue
            else:
                c_path_from = path_list[c_path_from_idx]
            event_indexed_segments = hist_entry[3]
            event_typed_segments = indexed_segments_to_typed_segments(event_indexed_segments, old_index_to_segment_dict)
            wt_segments, mt_segments, size_dict = form_mt_wt_indexed_paths_for_interpreter(c_path_from, new_segment_to_index_dict)
            interpreter_output, _ = interpret_haplotypes([mt_segments], [wt_segments], [c_path_from.path_chr], size_dict)
            current_deletion_info = ''
            for interpreted_event in interpreter_output:
                if interpreted_event[1] == 'deletion':
                    c_event_segments = interpreted_event[2][0].split('.')[2].replace('wt(', '').replace(')', '').split(',')
                    c_event_typed_segments = indexed_segments_to_typed_segments(c_event_segments, new_index_to_segment_dict)
                    if event_typed_segments == c_event_typed_segments:
                        current_deletion_info = interpreted_event[2][0]
                        break
            if current_deletion_info == '':
                # no deletion of the same segment found by interpreter
                sv_edges_list.append([])
                continue
            left_boundary_indexed_seg = current_deletion_info.split('.')[3]
            right_boundary_indexed_seg = current_deletion_info.split('.')[4]
            seg1 = indexed_segments_to_typed_segments([left_boundary_indexed_seg], new_index_to_segment_dict)[0]  # always only 1 seg
            seg4 = indexed_segments_to_typed_segments([right_boundary_indexed_seg], new_index_to_segment_dict)[0]
            c_sv_edges.append((seg1.chr_name, seg1.end, seg4.chr_name, seg4.start))
        elif event_type == 'nonreciprocal translocation':
            if (c_path_to_idx is not None) and (c_path_from_idx is not None):
                c_path_to = path_list[c_path_to_idx]
                c_path_from = path_list[c_path_from_idx]
                indexed_event_segments = hist_entry[3]
                event_typed_segments = indexed_segments_to_typed_segments(indexed_event_segments, old_index_to_segment_dict)

                # 1) resolve insertion in path_to
                sublist_idx = find_seg_sublist(c_path_to.linear_path.segments, event_typed_segments)
                if sublist_idx is not None:
                    seg1 = c_path_to.linear_path.segments[sublist_idx - 1]
                    seg2 = c_path_to.linear_path.segments[sublist_idx]
                    seg3 = c_path_to.linear_path.segments[sublist_idx + len(indexed_event_segments) - 1]
                    seg4 = c_path_to.linear_path.segments[sublist_idx + len(indexed_event_segments)]
                    c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                    c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
                else:
                    # make sure we don't write only half of the case, so we can manually distinguish the whole thing
                    sv_edges_list.append(c_sv_edges)
                    continue

                # 2) resolve deletion in path_from
                wt_segments, mt_segments, size_dict = form_mt_wt_indexed_paths_for_interpreter(c_path_from, new_segment_to_index_dict)
                interpreter_output, _ = interpret_haplotypes([mt_segments], [wt_segments], [c_path_from.path_chr], size_dict)
                current_deletion_info = ''
                for interpreted_event in interpreter_output:
                    if interpreted_event[1] == 'deletion':
                        c_event_segments = interpreted_event[2][0].split('.')[2].replace('wt(', '').replace(')', '').split(',')
                        c_event_typed_segments = indexed_segments_to_typed_segments(c_event_segments, new_index_to_segment_dict)
                        if event_typed_segments == c_event_typed_segments:
                            current_deletion_info = interpreted_event[2][0]
                            break
                if current_deletion_info != '':
                    left_boundary_indexed_seg = current_deletion_info.split('.')[3]
                    right_boundary_indexed_seg = current_deletion_info.split('.')[4]
                    seg1 = indexed_segments_to_typed_segments([left_boundary_indexed_seg], new_index_to_segment_dict)[0]  # always only 1 seg
                    seg4 = indexed_segments_to_typed_segments([right_boundary_indexed_seg], new_index_to_segment_dict)[0]
                    c_sv_edges.append((seg1.chr_name, seg1.end, seg4.chr_name, seg4.start))
                else:
                    # no deletion of the same segment found by interpreter
                    # make sure we don't write only half of the case, so we can manually distinguish the whole thing
                    sv_edges_list.append([])
                    continue
            else:
                sv_edges_list.append(c_sv_edges)
                continue
        elif event_type == 'arm tandem duplication':
            sv_edges_list.append(c_sv_edges)
            continue
        elif event_type == 'arm deletion':
            sv_edges_list.append(c_sv_edges)
            continue
        else:
            raise RuntimeError('event type undefined:', event_type)
        sv_edges_list.append(c_sv_edges)
    return sv_edges_list


def get_kt_index(input_index_dict, input_segment):
    for key in input_index_dict:
        if input_segment.same_segment_ignore_dir(key):
            return input_index_dict[key]


def get_t2_segments(masking_file):
    masking_arm = read_forbidden_regions(masking_file)
    t2_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == "telomere2":
            t2_segments.append(segment)
    return Arm(t2_segments, "t2_segments")


def segments_are_continuous(segment1: Segment, segment2: Segment):
    if segment1.direction():
        if segment1.end + 1 == segment2.start:
            return True
    else:
        if segment1.end - 1 == segment2.start:
            return True
    return False


def read_history_edges_intermediate_file(file_path, ignore_terminal_composite=True):
    event_sv_edges = []
    event_id = 0
    with open(file_path) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split(': ')
            history_label = line[0]
            edges = line[1]
            if edges.startswith('[\''):
                # chr_deletion OR chr_wide_event
                continue
            edges = eval(edges)
            history_label = eval(history_label)[0]
            if ignore_terminal_composite:
                history_label = history_label.replace('composite-', '').replace('-terminal', '')
            event_type = f"{history_label}#{event_id}"
            event_sv_edges.append({'label': event_type, 'edges': edges})
            event_id += 1
    return event_sv_edges


def history_edges_intermediate_file_terminal_event_labeling(input_dir_path, output_dir_path, forbidden_region_file):
    forbidden_df = pd.read_csv(forbidden_region_file, sep='\t')
    for file in os.listdir(input_dir_path):
        output_str = []
        with open(input_dir_path + file) as fp_read:
            for line in fp_read:
                line = line.replace('\n', '').split(': ')
                history_label = line[0]
                edges = line[1]
                edges = eval(edges)
                history_label = eval(history_label)
                if type(edges[0]) is str:
                    # chromosomal events
                    output_str.append('{}: {}'.format(history_label, edges))
                else:
                    event_is_terminal = False
                    for sv_edge in edges:
                        if edge_in_forbidden_region(str(sv_edge[0]),
                                                    int(sv_edge[1]),
                                                    str(sv_edge[2]),
                                                    int(sv_edge[3]),
                                                    forbidden_df):
                            event_is_terminal = True
                            break
                    if event_is_terminal:
                        new_event_type = history_label[0] + '-terminal'
                        history_label = (new_event_type, history_label[1], history_label[2], history_label[3])
                    output_str.append('{}: {}'.format(history_label, edges))
        output_str = '\n'.join(output_str)
        with open(output_dir_path + file, 'w') as fp_write:
            fp_write.write(output_str)


def edge_in_forbidden_region(chrom1, pos1, chrom2, pos2, forbidden_df):
    for index, row in forbidden_df.iterrows():
        if row['Chr'] == chrom1:
            if row['StartPos'] - 3 <= pos1 <= row['EndPos'] + 3:
                # pos1 is in forbidden region
                return True
        elif row['Chr'] == chrom2:
            if row['StartPos'] - 3 <= pos2 <= row['EndPos'] + 3:
                return True
    return False


def edge_conversion(indexed_seg1, indexed_seg2, index_to_segment_dict):
    typed_seg1 = index_to_segment_dict[indexed_seg1[:-1]]
    typed_seg2 = index_to_segment_dict[indexed_seg2[:-1]]
    chr1 = typed_seg1.chr_name
    chr2 = typed_seg2.chr_name
    if indexed_seg1[-1] == "-":
        pos1 = typed_seg1.start
    else:
        pos1 = typed_seg1.end
    if indexed_seg2[-1] == '-':
        pos2 = typed_seg2.end
    else:
        pos2 = typed_seg2.start
    print('(\'{}\', {}, \'{}\', {})'.format(chr1, pos1, chr2, pos2))


def edges_are_similar(edge1: (str, int, str, int), edge2: (str, int, str, int), allowance=5):
    if edge1[0] == edge2[0] and edge1[2] == edge2[2]:
        if abs(edge1[1] - edge2[1]) <= allowance and abs(edge1[3] - edge2[3]) <= allowance:
            return True
    return False

from KarUtils import *


def read_cluster_file(file_path):
    index_to_segment = {}
    karsim_path_list = []
    omkar_path_list = []
    with open(file_path) as fp_read:
        fp_read.readline()

        # fill segment dict
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            line = line.replace('\n', '').split('\t')
            index_to_segment[int(line[0])] = Segment(line[1], int(line[2]), int(line[3]), segment_type=line[4])

        # karsim paths
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            line = line.replace('\n', '').split('\t')
            segment_list = segment_indices_to_segments(line[2].split(','), index_to_segment)
            new_path = Path(Arm(segment_list, 'path'), path_name=line[0], path_chr=line[1])
            karsim_path_list.append(new_path)

        # omkar paths
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            line = line.replace('\n', '').split('\t')
            segment_list = segment_indices_to_segments(line[2].split(','), index_to_segment)
            new_path = Path(Arm(segment_list, 'path'), path_name=line[0], path_chr=line[1])
            omkar_path_list.append(new_path)

        # edges of interest
        edges_of_interest = []
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            line = line.replace('\n', '')
            info = line.split('\t')[0].replace('(', '').replace(')', '').split(', ')
            edge_distance = int(line.split('\t')[1])
            edges_of_interest.append((info[0], info[1], int(info[2]), info[3], edge_distance))

    return index_to_segment, karsim_path_list, omkar_path_list, edges_of_interest

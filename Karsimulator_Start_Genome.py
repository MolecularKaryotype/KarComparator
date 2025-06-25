import re

from KarUtils import *


def generate_genome_from_KT(input_file: str) -> Genome:
    """
    directly generating Genome from a formatted KT file
    :param input_file: path to input file
    :return: Genome Object
    """
    full_KT = {}
    motherboard_segments = []
    centromere_segments = []
    chromosome_of_interest = set()
    segment_dict = {}

    with open(input_file) as fp_read:
        line = fp_read.readline()  # column headers
        ## read-in segment indexing
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            info = line.replace('\n', '').split('\t')
            new_segment = Segment(info[1], int(info[2]), int(info[3]))
            segment_dict[info[0]] = new_segment
            chromosome_of_interest.add(info[1])
            if info[0][0] == "C":
                centromere_segments.append(new_segment)
            else:
                motherboard_segments.append(new_segment)

        line = fp_read.readline()  # column headers

        # construct slots
        for chromosome in chromosome_of_interest:
            full_KT[chromosome] = []

        ## read-in the segments in each chromosome
        while True:
            line = fp_read.readline().replace('\n', '')
            if line[0] == '-':
                break
            info = line.split('\t')

            if info[1] == 'deleted':
                new_chromosome = Chromosome(info[0], Arm([], 'deleted'), Arm([], 'deleted'), int(info[2]),
                                            int(info[3]), Arm([], 'deleted'), True)
                full_KT[info[0][:-1]].append(new_chromosome)
                continue

            p_arm_segments = []
            q_arm_segments = []
            current_centromere_segments = []
            segment_indices = info[1].split(',')
            # p_arm
            current_index = 0
            while segment_indices[current_index][0] != 'C':
                index_value = int(segment_indices[current_index][:-1])
                index_direction = segment_indices[current_index][-1]
                current_segment = segment_dict[str(index_value)].duplicate()
                if index_direction == "-":
                    current_segment.invert()
                p_arm_segments.append(current_segment)
                current_index += 1
            # centromere
            index_direction = segment_indices[current_index][-1]
            current_segment = segment_dict[segment_indices[current_index][:-1]].duplicate()
            if index_direction == "-":
                current_segment.invert()
            current_centromere_segments.append(current_segment)
            current_index += 1
            # q_arm
            while current_index < len(segment_indices):
                index_value = int(segment_indices[current_index][:-1])
                index_direction = segment_indices[current_index][-1]
                current_segment = segment_dict[str(index_value)].duplicate()
                if index_direction == "-":
                    current_segment.invert()
                q_arm_segments.append(current_segment)
                current_index += 1

            new_chromosome = Chromosome(info[0], Arm(p_arm_segments, 'p'), Arm(q_arm_segments, 'q'),
                                        int(info[2]), int(info[3]),
                                        Arm(current_centromere_segments, 'centromere'))
            # append chromosome to slot in genome, ignore the last char in chromosome's unique ID
            full_KT[info[0][:-1]].append(new_chromosome)

        ## read-in histories
        histories = []
        for line in fp_read:
            if line.startswith('block'):
                # first block now is reached
                break

        for line in fp_read:
            if line.startswith('block'):
                continue
            line = line.replace('\n', '')
            if len(line) == 0:
                # last line reached
                break
            event = line.split(', ')[0].split(' on')[0].replace('\t', '')
            event_segments = line.split('[')[1].split(']')[0].split(',')
            chr_info = line.split(', ')[1]
            chrs = re.findall(r'from\s+(.*?)\s+to\s+(.*?)$', chr_info)[0]
            chr_from = chrs[0]
            chr_to = chrs[1]
            current_history = (event, chr_from, chr_to, event_segments)
            histories.append(current_history)

    genome = Genome(full_KT, motherboard_segments, centromere_segments, histories)

    return genome


def get_history_events(input_file: str, origin_chr: [str]):
    # TODO: archive function, start_genome now has this functionality
    selected_chr_events = {}
    with open(input_file) as fp_read:
        section_number = 0
        for line in fp_read:
            if line.startswith('-'):
                section_number += 1
                if section_number == 2:
                    break
        for line in fp_read:
            if line.startswith('block'):
                break

        # first block now is reached
        for line in fp_read:
            if line.startswith('block'):
                continue

            line = line.replace('\n', '')
            if len(line) == 0:
                # last line reached
                break

            event = line.split(', ')[0].split(' on')[0].replace('\t', '')
            chr_info = line.split(', ')[1]
            chrs = re.findall(r'from\s+(.*?)\s+to\s+(.*?)$', chr_info)[0]
            chr_from = chrs[0][:-1]
            chr_to = chrs[1][:-1]
            if (chr_from in origin_chr) or (chr_to in origin_chr):
                if event in selected_chr_events:
                    selected_chr_events[event] += 1
                else:
                    selected_chr_events[event] = 1
    return selected_chr_events


def get_event_chr(input_file: str):
    # TODO: this function is redundant as we can now directly observe using the history
    event_chr = set()
    with open(input_file) as fp_read:
        section_number = 0
        for line in fp_read:
            if line.startswith('-'):
                section_number += 1
                if section_number == 2:
                    break
        for line in fp_read:
            if line.startswith('block'):
                break

        # first block now is reached
        for line in fp_read:
            if line.startswith('block'):
                continue

            line = line.replace('\n', '')
            if len(line) == 0:
                # last line reached
                break

            info = line.split(', ')[1]
            chrs = re.findall(r'from\s+(.*?)\s+to\s+(.*?)$', info)
            event_chr.add(chrs[0][0])
            event_chr.add(chrs[0][1])

    def custom_sort(input_chr):
        chr_index = input_chr[3:]
        chr_index = chr_index[:-1]
        if chr_index == "X":
            chr_index = 23
        elif chr_index == "Y":
            chr_index = 24
        else:
            chr_index = int(chr_index)

        return chr_index

    return sorted(list(event_chr), key=custom_sort)

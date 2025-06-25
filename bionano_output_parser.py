import pandas as pd


## SMAP
def smap_to_df(smap_filepath):
    with open(smap_filepath) as fp_read:
        line_count = 0
        for line in fp_read:
            line_count += 1
            if line.startswith('#h Smap'):
                break
        headers = line.replace('\n', '').split('\t')[1:]
    df = pd.read_csv(smap_filepath, sep='\t', skiprows=line_count+1, header=None, names=headers)
    return df


def filter_confidence(input_df, threshold=0.0):
    return input_df[(input_df['Confidence'] >= threshold) | (input_df['Type'].str.contains('translocation'))]


def filter_chroms(input_df, chrom1, chrom2):
    if chrom2 == -1 or chrom1 == chrom2:
        return input_df[(input_df['RefcontigID1'] == chrom1) | (input_df['RefcontigID2'] == chrom1)]
    else:
        return input_df[((input_df['RefcontigID1'] == chrom1) & (input_df['RefcontigID2'] == chrom2)) | ((input_df['RefcontigID2'] == chrom1) & (input_df['RefcontigID1'] == chrom2))]


def filter_event_type(input_df, event_types):
    return input_df[input_df['Type'].isin(event_types)]


def approx_edge_search(input_df, pos1, pos2, distance=20000, ignore_orientation=True):
    if ignore_orientation:
        mask1 = abs(input_df['RefStartPos'] - pos1) + abs(input_df['RefEndPos'] - pos2) < distance
        mask2 = abs(input_df['RefStartPos'] - pos2) + abs(input_df['RefEndPos'] - pos1) < distance
        return input_df[mask1 | mask2]


def approx_edge_search_inv(input_df, pos1, pos2, distance):
    filtered_df = input_df[input_df['Type'].str.contains('inversion', case=False)]
    if pos1 != pos2:
        mask1 = abs(filtered_df['RefStartPos'] - pos1) < distance
        mask2 = abs(filtered_df['RefEndPos'] - pos1) < distance
        mask3 = abs(filtered_df['RefStartPos'] - pos2) < distance
        mask4 = abs(filtered_df['RefEndPos'] - pos2) < distance
        return filtered_df[mask1 | mask2], filtered_df[mask3 | mask4]
    else:
        mask1 = abs(filtered_df['RefStartPos'] - pos1) < distance
        mask2 = abs(filtered_df['RefEndPos'] - pos1) < distance
        return filtered_df[mask1 | mask2], filtered_df[mask1 | mask2]

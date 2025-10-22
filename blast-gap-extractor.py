#------------------------------------------------------------------
# Blast-Gap-Extractor
#------------------------------------------------------------------
# 2025/06/18
# Motohiro Akashi
# motohiro-akashi@st.seikei.ac.jp
#------------------------------------------------------------------

import pandas as pd

# 1. loading BLASTn result data
blast_columns = [
    'query', 'subject', '% identity', 'alignment length', 'mismatches', 'gap opens',
    'q_start', 'q_end', 's_start', 's_end', 'evalue', 'bit score'
]
blast_df = pd.read_csv('all_blastn_out.tsv', sep='\t', names=blast_columns)

# 2. loading genome length data
genome_columns = ['acc.ver', 'length', 'gc%', 'gcsi']
genome_df = pd.read_csv('all_n_seqinfo.tsv', sep='\t', names=genome_columns)

# 3. remove self-comparisons (delete lines where query and subject are the same)
blast_df = blast_df[blast_df['query'] != blast_df['subject']]

# 4. dictionary of genome lengths
genome_lengths = dict(zip(genome_df['acc.ver'], genome_df['length']))

# 5. function to find undetected areas
def find_uncovered_regions(blast_df, genome_lengths):
    uncovered_regions = []

    # Integrate all coordinates to process query and subject in a unified manner
    for genome in genome_lengths.keys():
        length = genome_lengths[genome]

        # Get homology regions in both query and subject
        sub_df = blast_df[(blast_df['query'] == genome)][['q_start', 'q_end']].rename(columns={'q_start': 'start', 'q_end': 'end'})
        sub_df = pd.concat([sub_df,
                            blast_df[(blast_df['subject'] == genome)][['s_start', 's_end']].rename(columns={'s_start': 'start', 's_end': 'end'})])

        # If no data is available, the entire area is assumed to be uncovered.
        if sub_df.empty:
            uncovered_regions.append({'genome': genome, 'start': 1, 'end': length})
            continue

        # Sort by starting point and merge duplicates
        sub_df = sub_df.sort_values('start').reset_index(drop=True)
        merged = []
        current_start, current_end = sub_df.iloc[0]['start'], sub_df.iloc[0]['end']

        for _, row in sub_df.iloc[1:].iterrows():
            if row['start'] <= current_end + 1:  # Overlap or adjacent.
                current_end = max(current_end, row['end'])
            else:  # Non-contiguous regions.
                merged.append((current_start, current_end))
                current_start, current_end = row['start'], row['end']

        merged.append((current_start, current_end))  # Add last region.

        # Get uncovered regions.
        prev_end = 0
        for start, end in merged:
            if start > prev_end + 1:
                uncovered_regions.append({'genome': genome, 'start': prev_end + 1, 'end': start - 1})
            prev_end = end

        # Add if the end of the genome is uncovered.
        if prev_end < length:
            uncovered_regions.append({'genome': genome, 'start': prev_end + 1, 'end': length})

    # Rearrange uncovered_regions to DataFrame.
    df = pd.DataFrame(uncovered_regions)

    # Remove data if the whole genomic region is assigned as uncovered region
    if not df.empty:
        df = df[~(
            (df['start'] == 1) &
            (df.apply(lambda x: x['end'] == genome_lengths.get(x['genome'], -1), axis=1))
        )]

    return df

# 6. retrieve uncovered areas
uncovered_df = find_uncovered_regions(blast_df, genome_lengths)

# 7. save the results
uncovered_df.to_csv('uncovered_regions.tsv', sep='\t', index=False)

# 6. retrieve uncovered areas
uncovered_df = find_uncovered_regions(blast_df, genome_lengths)

# 7. save the results
if uncovered_df.empty:
    print("No uncovered regions detected in any genome.")
else:
    uncovered_df.to_csv('uncovered_regions.tsv', sep='\t', index=False)
    print("Uncovered regions have been saved to 'uncovered_regions.tsv'.")
    print(uncovered_df)

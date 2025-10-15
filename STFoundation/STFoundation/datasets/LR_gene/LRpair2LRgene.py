## 2025.10.15 convert LRpair to LRgene and save as csv file

#!/usr/bin/env python3
import os
import pandas as pd
import re
import argparse
import sys

def LRpair2LRgene(species='human', version='CellChatDBv2'):
    """
    Extract LR genes from CellChatDB interaction file and save as CSV.

    Args:
        species (str): Species name, default 'human'.
        version (str): Database version, default 'CellChatDBv2'.

    Returns:
        pd.DataFrame: Unique LR genes DataFrame.
    """
    dataset_dir = './STFoundation/datasets'
    if not os.path.exists(dataset_dir):
        print(f"Dataset directory not found: {dataset_dir}")
        sys.exit(1)
    os.chdir(dataset_dir)

    # Load interaction database
    interaction_file = f'LR_data/{species}-{version}_interaction.csv'
    if not os.path.exists(interaction_file):
        print(f"Interaction file not found: {interaction_file}")
        sys.exit(1)
    geneInter = pd.read_csv(interaction_file, header=0, index_col=0)
    print("CellChatDBv2 LRI dim:", geneInter.shape)
    print("CellChatDB_v1 Number:", sum(geneInter['version'] == 'CellChatDB v1'))
    print("CellChatDB_v2 Number:", sum(geneInter['version'] == 'CellChatDB v2'))

    # Display ligand and receptor arrays
    ligand = geneInter.ligand.values
    receptor = geneInter.receptor.values
    print(ligand)
    print(receptor)

    # Extract 'interaction_name_2' and split into genes
    geneInterPair = geneInter['interaction_name_2']
    split_genes = []
    for pair in geneInterPair:
        pair = pair.replace('(', '').replace(')', '')  # Remove parentheses
        genes = re.split(' - |:|\+', pair)             # Split by multiple delimiters
        split_genes.extend(genes)

    # Get unique genes and save to DataFrame
    LRgene = pd.Series(split_genes).unique()
    LRgene_df = pd.DataFrame(LRgene, columns=["LR gene"])
    print("LRgene_df shape:", LRgene_df.shape)
    print("CellChatDB_v1 Gene Number:", 963)
    print("CellChatDB_v2 Gene Number:", LRgene_df.shape[0])

    # Save as CSV (no row names)
    output_dir = "LR_gene"
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, f"LRgene_{version}_baseline_{species}.csv")
    LRgene_df.to_csv(output_file, index=False)
    print(f"Saved LR genes to: {output_file}")

    return LRgene_df

def main():
    parser = argparse.ArgumentParser(description="Extract LR genes from CellChatDB interaction file and save as CSV.")
    parser.add_argument('--species', type=str, default='human', choices=['human', 'mouse'],
                        help='Species name (default: human)')
    parser.add_argument('--version', type=str, default='CellChatDBv2', choices=['CellChatDBv2', 'CellChatDB'],
                        help='Database version (default: CellChatDBv2)')
    args = parser.parse_args()

    LRpair2LRgene(args.species, args.version)

if __name__ == '__main__':
    main()

## usage:
# python STFoundation/datasets/LR_gene/LRpair2LRgene.py --species mouse


import os
import argparse
import pandas as pd
import numpy as np
from ete3 import Tree
from Bio import SeqIO


def main(input_tree, input_fasta, input_tsv, output_tree, output_ali, output_traits):
    assert (os.path.isfile(input_tree))
    assert (os.path.isfile(input_fasta))
    assert (os.path.isfile(input_tsv))

    os.makedirs(os.path.dirname(output_ali), exist_ok=True)
    tree = Tree(input_tree)
    fasta_dico = {f.id: f.seq for f in SeqIO.parse(open(input_fasta, 'r'), 'fasta')}
    print(f"Found {len(fasta_dico)} sequences in {input_fasta}.")
    if len(fasta_dico) == 0:
        exit(0)
    assert set(fasta_dico) == set(tree.get_leaf_names()), f"{list(fasta_dico)} \n {list(tree.get_leaf_names())}"

    df_traits = pd.read_csv(input_tsv, sep="\t")
    print(f"Found {len(df_traits)} species in {input_tsv}.")

    assert set(fasta_dico) == set(df_traits["TaxonName"]), f"{list(fasta_dico)} \n {list(df_traits['TaxonName'])}"
    species_to_keep = set(tree.get_leaf_names()).intersection(set(df_traits["TaxonName"])).intersection(set(fasta_dico))
    print(f"Found {len(species_to_keep)} species in common between the tree, fasta and traits files.")

    # Write tree file
    tree.prune(species_to_keep, preserve_branch_length=True)
    assert (len(tree.get_leaf_names()) == len(species_to_keep))
    tree.write(outfile=output_tree, format=1)

    # Write ali file
    ids_seqs = [(f_id, str(fasta).replace("!", "-").replace("?", "-")) for f_id, fasta in fasta_dico.items() if
                f_id in species_to_keep]
    assert len(ids_seqs) == len(tree.get_leaf_names()), f"{len(ids_seqs)} != {len(tree.get_leaf_names())}"
    ali_file = open(output_ali, 'w')
    ali_file.write(str(len(ids_seqs)) + " " + str(len(ids_seqs[0][1])) + "\n")
    ali_file.write("\n".join([" ".join(id_seq) for id_seq in ids_seqs]))
    ali_file.close()

    # Write traits file
    df_traits.to_csv(output_traits, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--input_tree", type=str, required=True, help="Input tree file")
    parser.add_argument("-f", "--input_fasta", type=str, required=True, help="Input fasta file")
    parser.add_argument("-p", "--input_tsv", type=str, required=True, help="Input tsv file")
    parser.add_argument("-o", "--output_tree", type=str, required=True, help="Output tree file")
    parser.add_argument("-a", "--output_ali", type=str, required=True, help="Output alignment file")
    parser.add_argument("-r", "--output_traits", type=str, required=True, help="Output traits file")
    args = parser.parse_args()
    main(args.input_tree, args.input_fasta, args.input_tsv, args.output_tree, args.output_ali, args.output_traits)


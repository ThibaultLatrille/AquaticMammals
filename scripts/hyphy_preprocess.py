import os
import argparse
import pandas as pd
from ete3 import Tree
from Bio import SeqIO


def write_fasta(dico_fasta, output, species=None):
    outfile = open(output, "w")
    outfile.write("\n".join([f">{seq_id}\n{seq}" for seq_id, seq in dico_fasta.items() if
                             (species is None) or (seq_id in species)]))
    outfile.close()


def main(input_tree, input_fasta, input_tsv, output_tree, output_ali):
    assert (os.path.isfile(input_tree))
    assert (os.path.isfile(input_fasta))
    assert (os.path.isfile(input_tsv))

    os.makedirs(os.path.dirname(output_ali), exist_ok=True)
    tree = Tree(input_tree)
    print(f"Found {len(tree.get_leaf_names())} species in {input_tree}.")
    fasta_dico = {f.id: f.seq for f in SeqIO.parse(open(input_fasta, 'r'), 'fasta')}
    print(f"Found {len(fasta_dico)} sequences in {input_fasta}.")
    if set(fasta_dico) != set(tree.get_leaf_names()):
        set_diff = set(fasta_dico).symmetric_difference(set(tree.get_leaf_names()))
        print(f"Warning: {set_diff} are different between the fasta and the tree files.")

    df_traits = pd.read_csv(input_tsv, sep="\t")
    print(f"Found {len(df_traits)} species in {input_tsv}.")
    if set(fasta_dico) != set(df_traits["TaxonName"]):
        set_diff = set(fasta_dico).symmetric_difference(set(df_traits["TaxonName"]))
        print(f"Warning: {set_diff} are different between the fasta and the traits files.")

    species_to_keep = set(tree.get_leaf_names()).intersection(set(df_traits["TaxonName"])).intersection(set(fasta_dico))
    print(f"Found {len(species_to_keep)} species in common between the tree, fasta and traits files.")
    if len(species_to_keep) == 0:
        exit(1)

    # Write tree file
    tree.prune(species_to_keep, preserve_branch_length=True)
    assert (len(tree.get_leaf_names()) == len(species_to_keep))
    # Find all the leaves in the tree that are aquatic species
    aquatic_species = df_traits[df_traits["Aquatic_adaptation"] == 4]["TaxonName"].values

    # Set to T if all the leaves of the subtree are aquatic species, R otherwise
    for node in tree.traverse('preorder'):
        if all([leaf.name in aquatic_species for leaf in node.get_leaves()]):
            node.name  = node.name + "{T}"
        else:
            node.name  = node.name + "{R}"

    print(f"Found {len([leaf for leaf in tree.get_leaves() if "{T}" in leaf.name])} aquatic leaves in the tree.")
    print(f"Found {len([node for node in tree.traverse() if "{T}" in node.name])} aquatic nodes in the tree.")
    if len([leaf for leaf in tree.get_leaves() if "{R}" in leaf.name]) == 0:
        exit(1)
    tree.write(outfile=output_tree, format=1)

    # Write fasta file
    write_fasta(fasta_dico, output_ali, species=species_to_keep)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--input_tree", type=str, required=True, help="Input tree file")
    parser.add_argument("-f", "--input_fasta", type=str, required=True, help="Input fasta file")
    parser.add_argument("-p", "--input_tsv", type=str, required=True, help="Input tsv file")
    parser.add_argument("-o", "--output_tree", type=str, required=True, help="Output tree file")
    parser.add_argument("-a", "--output_ali", type=str, required=True, help="Output alignment file")
    args = parser.parse_args()
    main(args.input_tree, args.input_fasta, args.input_tsv, args.output_tree, args.output_ali)

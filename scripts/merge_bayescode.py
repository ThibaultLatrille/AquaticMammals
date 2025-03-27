#!/usr/bin/env python3
import os
from collections import defaultdict
from ete3 import Tree
import argparse
import pandas as pd
import numpy as np

dico_coeff = {
    "covariances": "cov",
    "correlation coefficients": "cor",
    "posterior probabilities of a positive coefficient": "ppos",
    "precisions": "prec",
    "partial correlation coefficients": "pcor",
    "posterior probabilities of a positive partial coefficient": "ppos_pcor"
}


def extract_matrix(header, f, name):
    matrix = [[0. for _ in range(len(header))] for _ in range(len(header))]
    assert f.readline().strip() == name
    assert f.readline().strip() == ""
    row = 0
    for line in f:
        if line.strip() == "":
            break
        row_list = [x.strip() for x in line.strip().split(" ") if x.strip() != ""]
        matrix[row] = [float(x) if x not in ["-"] else x for x in row_list]
        assert len(matrix[row]) == len(header)
        row += 1
    return matrix


def open_covar_file(covar_file):
    """
    The file has the following format:
    entries are in the following order:
    Phenotype_mean
    Genotype_mean

    covariances

    62.266  61.58
    61.58 61.114

    correlation coefficients

    1  0.998
    0.998      1

    posterior probs

    -      1
    1      -
    """
    matrices = {}
    header = []
    with open(covar_file) as f:
        f.readline()
        for line in f:
            if line.strip() == "":
                break
            header.append(line.strip())
        for m in dico_coeff.keys():
            matrices[m] = extract_matrix(header, f, m)
    return header, matrices


def main(input_data, output_cov, output_omega, input_traits):
    omega_dict = defaultdict(lambda: defaultdict(lambda: "NA"))
    output_cov_dict = defaultdict(list)

    species_list = set()
    cds_list = sorted(os.listdir(input_data))
    for cds in cds_list:
        cds_path = f"{input_data}/{cds}"
        if not os.path.isdir(cds_path) or cds.startswith("."):
            continue

        omega_path = f"{cds_path}/nodeomega_1.Omega.nhx"
        if not os.path.exists(omega_path):
            print(f"Missing {omega_path}")
            continue

        omega_tree = Tree(omega_path, format=1)
        for leave in omega_tree.get_leaves():
            species_list.add(leave.name)
            omega_dict[cds][leave.name] = leave.Omega

        omega_list = [float(node.Omega) for node in omega_tree.traverse() if "Omega" in node.features]
        nb_omega_geq1 = sum([1 for omega in omega_list if omega >= 1])
        mean_omega = np.mean(omega_list)

        cov_path = f"{cds_path}/nodeomega_1.cov"
        if not os.path.exists(cov_path):
            print(f"Missing {cov_path}")
            continue

        ali_path = f"{cds_path}/placnr.ali"
        nb_species, nb_sites = open(ali_path).readline().strip().split(" ")
        if input_traits != "":
            assert os.path.exists(input_traits)
            traits_path = input_traits
        else:
            traits_path = f"{cds_path}/placnr.traits"
        df_traits = pd.read_csv(traits_path, sep="\t")
        output_cov_dict["id"].append(cds)
        output_cov_dict["nb_species"].append(nb_species)
        output_cov_dict["nb_sites"].append(nb_sites)
        output_cov_dict["nb_branches"].append(len(omega_list))
        output_cov_dict["nb_omega_geq1"].append(nb_omega_geq1)
        output_cov_dict["mean_omega"].append(mean_omega)

        if "historical_Ne" in df_traits.columns:
            output_cov_dict["historical_Ne"].append(df_traits["historical_Ne"].mean())
            nb_ne = sum(np.isfinite(df_traits["historical_Ne"]))
            output_cov_dict["nb_ne"].append(nb_ne)

        header, matrices = open_covar_file(cov_path)
        for i in range(len(header)):
            for j in range(i, len(header)):
                if i == j:
                    continue
                for k, matrix in matrices.items():
                    output_cov_dict[f"{header[i]}_{header[j]}_{dico_coeff[k]}"].append(matrix[i][j])

    # Write cov file
    df_cov = pd.DataFrame(output_cov_dict)
    df_cov.to_csv(output_cov, sep="\t", index=False)

    # Write omega file
    sorted_species_list = sorted(species_list)
    output_omega_dict = defaultdict(list)
    for cds, omega in omega_dict.items():
        output_omega_dict["id"].append(cds)
        for species in sorted_species_list:
            output_omega_dict[species].append(omega[species])
    df_omega = pd.DataFrame(output_omega_dict)
    df_omega.to_csv(output_omega, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_traits", type=str, required=False, default="", help="Input traits file")
    parser.add_argument("--input_data", type=str, required=True, help="Folder containing data files")
    parser.add_argument("--output_cov", type=str, required=True, help="Output cov file")
    parser.add_argument("--output_omega", type=str, required=True, help="Output omega file")
    args = parser.parse_args()
    main(args.input_data, args.output_cov, args.output_omega, args.input_traits)

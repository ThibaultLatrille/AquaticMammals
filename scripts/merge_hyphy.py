#!/usr/bin/env python3
import os
from collections import defaultdict
import argparse
import json
import polars as pl


def main(input_data, output_file):
    results_dico = defaultdict(list)

    cds_list = sorted(os.listdir(input_data))
    for cds in cds_list:
        cds_path = f"{input_data}/{cds}"
        if not os.path.isdir(cds_path) or cds.startswith("."):
            continue
        json_path = f"{cds_path}/placnr.fasta.RELAX.json"
        if not os.path.exists(json_path):
            print(f"Missing {json_path}")
            continue
        # Read the JSON file with package json
        with open(json_path, 'r') as f:
            data = json.load(f)
            test = data["test results"]
            # Extract the relevant information
            results_dico["id"].append(cds)
            results_dico["LRT"].append(test["LRT"])
            results_dico["p-value"].append(test["p-value"])
            results_dico["k"].append(test["relaxation or intensification parameter"])

    # Save with polars
    pl_df = pl.DataFrame(results_dico, strict=False)
    # Sort the dataframe by the first column
    pl_df = pl_df.sort("id")
    pl_df.write_csv(output_file, separator="\t", include_header=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_data", type=str, required=True, help="Folder containing data files")
    parser.add_argument("--output_file", type=str, required=True, help="Output file")
    args = parser.parse_args()
    main(args.input_data, args.output_file)

#!/usr/bin/env python3
import os
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols


def match_refseq_to_ensembl_df(df_input, orthomam_summary):
    df_orthomam_summary = pd.read_csv(orthomam_summary, sep="\t")
    id_to_name = dict(zip(df_orthomam_summary["id"], df_orthomam_summary["name"]))
    id_to_ensg = dict(zip(df_orthomam_summary["id"], df_orthomam_summary["ensg"]))
    df_input["ensg"] = df_input["id"].map(id_to_ensg)
    df_input["name"] = df_input["id"].map(id_to_name)
    return df_input


def match_refseq_list_to_ensembl(refseq_list, orthomam_summary):
    df_orthomam_summary = pd.read_csv(orthomam_summary, sep="\t")
    id_to_ensg = dict(zip(df_orthomam_summary["id"], df_orthomam_summary["ensg"]))
    return [id_to_ensg[refseq] for refseq in refseq_list]


def mask_to_cds_set(input_mask):
    return set(pd.read_csv(input_mask, sep="\t")["id"])


def parse_df_omega(input_omega, orthomam_summary, species, input_mask):
    df_omega = pd.read_csv(input_omega, sep="\t", usecols=["id", species])
    mask_set = set(pd.read_csv(input_mask, sep="\t")["id"])
    print(f"Number of genes in OrthoMam: {len(df_omega)}")
    df_omega = df_omega[df_omega["id"].isin(mask_set)]
    print(f"Number of genes in OrthoMam with mask: {len(df_omega)}")
    df_omega = match_refseq_to_ensembl_df(df_omega, orthomam_summary)
    df_omega = df_omega.rename(columns={species: "omega"})
    # filter out genes with no omega
    df_omega = df_omega[df_omega["omega"].notnull()]
    print(f"Number of genes in OrthoMam with ENSG id: {len(df_omega)}")
    return df_omega


def merge_condition(condition):
    condition = condition.lower().strip()
    brain_parts = ["cortex", "hippocampus", "gyrus", "cerebellum", "thalamus", "lobe", "pole", "brain",
                   "caudate nucleus", "amygdala", "striatum", "substantia nigra", "pituitary gland", "pineal gland",
                   "corpus callosum"]
    if any([part in condition for part in brain_parts]):
        return "brain"
    heart_parts = ["atrium", "atrial", "ventricle", "heart", "mitral valve", "aorta", "aortic", "valve",
                   "pulmonary vein", "pulmonary artery", "ventricular septum", "vena cava"]
    if any([part in condition for part in heart_parts]):
        return "heart"
    colon_parts = ["colon", "rectum", "sigmoid", "cecum", "appendix", "anal canal", "anus", "rectal"]
    if any([part in condition for part in colon_parts]):
        return "colon"
    for organ in ["liver", "kidney", "testis"]:
        if organ in condition:
            return organ
    return "other"


def tf_log2TPM(df_array):
    k = np.log(2)
    return df_array.map(lambda x: x * k if np.isfinite(x) else x)


def tf_pLevel(df_array):
    sum_pLevel = np.nansum(df_array)
    return df_array.map(lambda x: np.log(x / sum_pLevel) if x > 0 else np.nan)


def parse_df_pLevel(input_pLevel, input_samples, species):
    df_pLevel = pd.read_csv(input_pLevel, sep="\t", dtype={"Gene ID": str, "Gene Name": str})
    df_pLevel.rename(columns={"Gene ID": "ensg", "Gene Name": "name"}, inplace=True)
    df_samples = pd.read_csv(input_samples, sep="\t")
    df_samples = df_samples[df_samples["organism"] == species]
    df_samples = df_samples[df_samples["disease"] == "normal"]
    df_samples["organism part"] = df_samples["organism part"].apply(lambda x: merge_condition(x))
    discard_conditions = {"other", "heart", "kidney"}
    df_samples = df_samples[~df_samples["organism part"].isin(discard_conditions)]
    df_samples = df_samples[df_samples["group"].isin(set(df_pLevel.columns))]
    # Remove organs with less than 3 samples
    df_samples = df_samples.groupby("organism part").filter(lambda x: len(x) > 5)

    set_conditions = set(df_samples["organism part"])
    sample_condition = {k: v for k, v in zip(df_samples["group"], df_samples["organism part"])}
    print(f"{len(set_conditions)} conditions: {set_conditions}")
    sample_set = sorted(set(df_samples["group"]).intersection(set(df_pLevel.columns)))
    return df_pLevel, sample_set, sample_condition


def regroup_pLevel(df_pLevel, sample_set, sample_condition):
    condition_to_sample = defaultdict(list)
    for sample in sample_set:
        condition_to_sample[sample_condition[sample]].append(sample)
    dico_output = dict()
    for condition, samples in condition_to_sample.items():
        dico_output[condition] = df_pLevel[["ensg"] + samples]
    return dico_output


def group_df_pLevel(input_pLevel, input_samples, species):
    df_pLevel, sample_set, sample_condition = parse_df_pLevel(input_pLevel, input_samples, species)
    return regroup_pLevel(df_pLevel, sample_set, sample_condition)


def group_df_log2TPM(input_expression_list):
    dico_output = dict()
    for input_expression in input_expression_list:
        condition = os.path.basename(os.path.dirname(input_expression)).split("_")[0]
        df_expression = pd.read_csv(input_expression, sep="\t")
        print(f"Log2TPM condition: {condition} with {len(df_expression)} CDS and {len(df_expression.columns)} samples")
        df_expression["ensg"] = df_expression.index
        dico_output[condition] = df_expression
    return dico_output


def compute_corr_stats(output_dict, df_ols, x, y):
    r_spearman = df_ols[x].corr(df_ols[y], method="spearman")
    output_dict["r_spearman"].append(r_spearman)
    r_pearson = df_ols[x].corr(df_ols[y], method="pearson")
    output_dict["r_pearson"].append(r_pearson)

    ols_res = ols(f"{y} ~ {x}", data=df_ols).fit()
    output_dict["slope"].append(ols_res.params[x])
    output_dict["slope_se"].append(ols_res.bse[x])
    output_dict["intercept"].append(ols_res.params["Intercept"])
    output_dict["intercept_se"].append(ols_res.bse["Intercept"])
    output_dict["r2"].append(ols_res.rsquared)
    output_dict["r2_adj"].append(ols_res.rsquared_adj)
    return ols_res


def y_from_x(x_val, ols_res, x):
    return ols_res.params["Intercept"] + ols_res.params[x] * x_val


def plot_scatter(df, x, y, ax, ols_res, title=None, x_label=None, y_label=None):
    ax.scatter(df[x], df[y], alpha=0.2)
    x_min, x_max = df[x].min(), df[x].max()
    ax.plot([x_min, x_max], [y_from_x(x_min, ols_res, x), y_from_x(x_max, ols_res, x)], color="red")
    if x_label is None:
        x_label = x.replace("_", " ")
    if y_label is None:
        y_label = y.replace("_", " ")
    if title is None:
        title = f"{y_label} ~ {x_label}\nR2: {ols_res.rsquared:.2f}, P-value: {ols_res.f_pvalue:.2e}"
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)


def plot_slopes(df, x, y, output, ols_res=None):
    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    if ols_res is None:
        ols_res = ols(f"{y} ~ {x}", data=df).fit()
    ax = axes[0]
    plot_scatter(df, x, y, ax, ols_res)
    ax = axes[1]
    # Bin the x values
    df["bin"] = pd.qcut(df[x], q=15, duplicates="drop")
    df_grouped = df.groupby("bin", observed=False)
    groupex = df_grouped[x].mean()
    groupey = df_grouped[y].mean()
    assert len(groupex) == len(groupey) == len(df_grouped)
    ax.scatter(groupex, groupey)
    x_min, x_max = groupex.min(), groupex.max()
    ax.plot([x_min, x_max], [y_from_x(x_min, ols_res, x), y_from_x(x_max, ols_res, x)], color="red")
    ax.set_title(f'Binned {x.replace("_", " ")} ~ {y.replace("_", " ")}')
    ax.set_xlabel(x.replace("_", " "))
    ax.set_ylabel(y.replace("_", " "))
    plt.tight_layout()
    plt.savefig(output)
    plt.close("all")
    plt.clf()

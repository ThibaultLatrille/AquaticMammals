import os

configfile: 'config/config.yaml'

FOLDER = os.path.abspath('.')
DATA_FOLDER = f"{FOLDER}/data/selected_genes"
PROC_DATA = f"{FOLDER}/data_processed"

CDS_list = set([i.split("_")[0] for i in os.listdir(DATA_FOLDER) if not i.startswith(".") and i.endswith(".fasta")])
# I modified "51725_NT_AL.tsv" to include the 191 species (file was malformed)
# I modified "55435_NT_AL.rootree" that I copied from OrthoMam v12
# For two CDS, the species don't match between files, I used the intersection: "3620_NT_AL"
print(f"{len(CDS_list)} CDS found in {DATA_FOLDER}.")
CDS_list = list(sorted(CDS_list))
if "CDS_LIST" in config:
    CDS_list = [i for i in CDS_list if i in map(str,config['CDS_LIST'])]
else:
    CDS_list = CDS_list[config['CDS_START_LIST']:config['CDS_END_LIST']]

hyphy_path = "/opt/homebrew/Caskroom/miniforge/base/envs/bioinfo/bin/hyphy"
if not os.path.exists(hyphy_path):
    # Find executable in the path using whereis. If not found, raise an error.
    split = os.popen(f'whereis hyphy').read().split()
    if len(split) > 1:
        hyphy_path = split[1].strip()
    else:
        raise FileNotFoundError(f'hyphy not found. Please install Hyphy and add it to your path.')

localrules: all,hyphy_preprocess,all_preprocess,merge_hyphy

rule all:
    input:
        f"{PROC_DATA}/Hyphy_merged.tsv"

rule hyphy_preprocess:
    input:
        scripts=f"{FOLDER}/scripts/hyphy_preprocess.py",
        fasta=f"{DATA_FOLDER}/{{CDS}}_NT_AL.fasta",
        tree=f"{DATA_FOLDER}/{{CDS}}_NT_AL.rootree",
        tsv=f"{DATA_FOLDER}/{{CDS}}_NT_AL.tsv"
    output:
        ali=f"{PROC_DATA}/Hyphy/{{CDS}}/placnr.fasta",
        tree=f"{PROC_DATA}/Hyphy/{{CDS}}/placnr.tree"
    shell:
        "python3 {input.scripts} --input_fasta {input.fasta} --input_tree {input.tree} --input_tsv {input.tsv} --output_ali {output.ali} --output_tree {output.tree}"

rule all_preprocess:
    input:
        expand(rules.hyphy_preprocess.output.tree,CDS=CDS_list)

rule run_hyphy:
    input:
        hyphy=hyphy_path,
        batch=f"{FOLDER}/hyphy-analyses/RELAX-mod/RELAX.bf",
        ali=rules.hyphy_preprocess.output.ali,
        tree=rules.hyphy_preprocess.output.tree
    output:
        log=f"{PROC_DATA}/Hyphy/{{CDS}}/hyphy.log",
        json=f"{PROC_DATA}/Hyphy/{{CDS}}/placnr.fasta.RELAX.json",
    params:
        time="2-23:00",mem=1000,name=lambda wildcards: f"run_CDS{wildcards.CDS}",
    shell: "{input.hyphy} {input.batch} --alignment {input.ali} --tree {input.tree} --test T --models Minimal > {output.log} 2>&1;"

rule merge_hyphy:
    input:
        script=f"{FOLDER}/scripts/merge_hyphy.py",
        json=expand(rules.run_hyphy.output.json,CDS=CDS_list)
    output:
        tsv=f"{PROC_DATA}/Hyphy_merged.tsv"
    shell:
        "python3 {input.script} --input_data {PROC_DATA}/Hyphy --output_file {output.tsv}"

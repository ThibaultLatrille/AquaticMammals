import os
import gzip


def open_file(path, rw="r"):
    return gzip.open(path, f'{rw}t') if path.endswith(".gz") else open(path, rw)


def open_fasta(path, filter_species=None) -> dict:
    if filter_species is None:
        filter_species = []
    outfile = {}
    ali_file = open_file(path)
    for seq_id in ali_file:
        sp_id = seq_id.replace('>', '').strip()
        seq = ali_file.readline().strip()
        if (len(filter_species) > 0) and (sp_id not in filter_species):
            continue
        outfile[sp_id] = seq
    return outfile


def write_fasta(dico_fasta, output):
    outfile = open_file(output, "w")
    outfile.write("\n".join([f">{seq_id}\n{seq}" for seq_id, seq in dico_fasta.items()]))
    outfile.close()

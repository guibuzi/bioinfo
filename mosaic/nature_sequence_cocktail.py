import random

def read_fasta(path):
    seq_dict = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                label = line[1: -1]
                seq_dict[label] = ''
            else:
                seq = line.replace("\n", "").strip()
                seq_dict[label] += seq
    return seq_dict


seqs = read_fasta("/home/zeng/python_work/bioinfo/mosaic/us_ha_sampling.fasta")

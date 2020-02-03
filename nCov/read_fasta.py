import os, sys
from collections import defaultdict


def read_fasta(_in):
    info = defaultdict(dict)
    with open(_in) as f:
        for line in f:
            if line.startswith(">"):
                attrs = line[1:-1].split("|")
                genbank_ac = attrs[0]
                info[genbank_ac]['strain_name'] = attrs[1]
                info[genbank_ac]['segment'] = attrs[2]
                info[genbank_ac]['date'] = attrs[3]
                info[genbank_ac]['host'] = attrs[4]
                info[genbank_ac]['country'] = attrs[5]
                info[genbank_ac]['subtype'] = attrs[6]
                info[genbank_ac]['virus_species'] = attrs[7]
                info[genbank_ac]['seq'] = ''
            else:
                info[genbank_ac]['seq'] += line.replace("\n", "")
    return info


def read_info(_in):
    info = defaultdict(dict)
    with open(_in) as f:
        for line in f:
            if line.startswith(">"):
                attrs = line[1:-1].split("|")
                genbank_ac = attrs[0]
                info[genbank_ac]['strain_name'] = attrs[1]
                info[genbank_ac]['segment'] = attrs[2]
                info[genbank_ac]['date'] = attrs[3]
                info[genbank_ac]['host'] = attrs[4]
                info[genbank_ac]['country'] = attrs[5]
                info[genbank_ac]['subtype'] = attrs[6]
                info[genbank_ac]['virus_species'] = attrs[7]
    return info


_229e = read_fasta("/home/zeng/Desktop/COV/SARS.fasta")
print(_229e)
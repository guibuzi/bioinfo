import random
import re


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


def find_crossover_point(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    strings_of_seq1 = set([seq1[i:i+8] for i in range(len1-8)])
    strings_of_seq2 = set([seq2[i:i+8] for i in range(len2-8)])
    common_strings = strings_of_seq1.intersection(strings_of_seq2)

    crossover1 = common_strings.pop()
    crossover2 = common_strings.pop()
    m1_1 = re.finditer(crossover1, seq1)
    m1_2 = re.finditer(crossover2, seq1)
    m2_1 = re.finditer(crossover1, seq2)
    m2_2 = re.finditer(crossover2, seq2)
    locations_crossover1_1 = [m1_i.start() for m1_i in m1_1]
    locations_crossover1_2 = [m2_i.start() for m2_i in m1_2]
    locations_crossover2_1 = [m1_i.start() for m1_i in m2_1]
    locations_crossover2_2 = [m2_i.start() for m2_i in m2_2]
    index_1_1 = random.choice(locations_crossover1_1)
    index_1_2 = random.choice(locations_crossover1_2)
    index_2_1 = random.choice(locations_crossover2_1)
    index_2_2 = random.choice(locations_crossover2_2)
    offset = random.randint(0,7)
    return (index_1_1+offset, index_1_2+offset, index_2_1+offset, index_2_2+offset)


def recombine(seq1, seq2, location):
    index1_1, index1_2, index2_1, index2_2 = location
    seq1_1 = seq1[:index1_1]
    seq1_2 = seq1[index1_1:index1_2]
    seq1_3 = seq1[index1_2:]
    seq2_1 = seq2[:index2_1]
    seq2_2 = seq2[index2_1:index2_2]
    seq2_3 = seq2[index2_2:]
    return (seq1_1+seq2_2+seq1_3, seq2_1+seq1_2+seq2_3)


def main(seq_list, k, n):
    populations = []
    for i in range(k):
        population_i = []
        for j in range(n):
            seq1 = random.choice(seq_list)
            seq2 = random.choice(seq_list)
            locations = find_crossover_point(seq1, seq2)
            child1, child2 = recombine(seq1, seq2, locations)
            population_i.append(child1)
            population_i.append(child2)
        populations.append(population_i)
    return populations


if __name__ == "__main__":
    data = read_fasta("mosaic/us_ha_sampling.fasta")
    seq_list = list(data.values())
    popu = main(seq_list, 4, len(seq_list))

import re, os, time
from collections import Counter
from multiprocessing import Process, Pool, Queue, Manager
import matplotlib.pyplot as plt


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


def score_single_location(epitope_list, k):
    total = len(epitope_list)
    count_by_epitope = Counter(epitope_list)
    top_k = count_by_epitope.most_common(k)
    coverage = [i[1] for i in top_k]
    return sum(coverage) / total


def score_k(to_score, k):
    length = len(to_score[0])
    score_list = []
    for i in range(length - 9):
        epitope_list = [seq[i:i+9] for seq in to_score]
        score = score_single_location(epitope_list, k)
        score_list.append(score)
    return score_list


def main():
    data = read_fasta("mosaic/us_ha_sampling_align.fasta")
    seq_list = list(data.values())
    scores_list = []
    for k in range(1, 9):
        scores = score_k(seq_list, k)
        scores_list.append((k, scores))
    return scores_list


if __name__ == "__main__":
    scores_list = main()
    for k in range(0, 8):
        plt.plot(scores_list[k][1], label=scores_list[k][0])
    plt.legend(loc=4)
    plt.title("The coverage of HA in U.S. 2009-2019")
    plt.savefig('test.jpg')

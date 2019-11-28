import os
import re
import time
from collections import Counter
from multiprocessing import Manager, Pool, Process, Queue

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


def main(path):
    data = read_fasta(path)
    seq_list = list(data.values())
    scores_list = []
    for k in range(1, 9):
        scores = score_k(seq_list, k)
        scores_list.append((k, scores))
    return scores_list


if __name__ == "__main__":
    scores_list_ha = main("mosaic/us_ha_sampling_align.fasta")
    scores_list_na = main("mosaic/us_na_sampling_align.fasta")
    scores_list_ha_sorted = [(x[0], sorted(x[1], reverse=True)) for x in scores_list_ha]
    scores_list_na_sorted = [(x[0], sorted(x[1], reverse=True)) for x in scores_list_na]

    fig = plt.figure(figsize=(10, 8))

    ax1 = fig.add_subplot(2, 2, 1)
    for k in range(0, 8):
        ax1.plot(scores_list_ha[k][1], label=scores_list_ha[k][0])
    plt.title("The coverage of HA in U.S. 2009-2019")  
    
    ax2 = fig.add_subplot(2, 2, 2)
    for k in range(0, 8):
        ax2.plot(scores_list_na[k][1], label=scores_list_na[k][0])
    plt.title("The coverage of NA in U.S. 2009-2019")

#    plt.legend(bbox_to_anchor=(1.2, 0.5), loc='right', fontsize='small')
#    fig.savefig('test.jpg')
 
#    fig = plt.figure(figsize=(10, 4))

    ax3 = fig.add_subplot(2, 2, 3)
    for k in range(0, 8):
        ax3.plot(scores_list_ha_sorted[k][1], label=scores_list_ha_sorted[k][0])
#    plt.title("The coverage of HA in U.S. 2009-2019")  
    
    ax4 = fig.add_subplot(2, 2, 4)
    for k in range(0, 8):
        ax4.plot(scores_list_na_sorted[k][1], label=scores_list_na_sorted[k][0])
#    plt.title("The coverage of NA in U.S. 2009-2019")

    plt.legend(bbox_to_anchor=(1.2, 0.5), loc='right', fontsize='small')
    fig.savefig('test2.jpg')

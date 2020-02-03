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
    length = len(to_score[1])
    score_list = []
    for i in range(length - 9):
        epitope_list = [seq[i:i+9] for seq in to_score]
        score = score_single_location(epitope_list, k)
        score_list.append(score)
    return score_list


def main(path):
    data = read_fasta(path)
    seq_list = list(set(list(data.values())))
    scores_list = []
    for k in range(1, 9):
        scores = score_k(seq_list, k)
        scores_list.append((k, scores))
    return scores_list


if __name__ == "__main__":
    score_list_h1 = main("/home/zeng/Desktop/ncbi_h1.fa")
    score_list_h3 = main("/home/zeng/Desktop/ncbi_h3.fa")
    score_list_n1 = main("/home/zeng/Desktop/ncbi_n1.fa")
    score_list_n2 = main("/home/zeng/Desktop/ncbi_n2.fa")

    score_list_h1_sorted = [(x[0], sorted(x[1], reverse=True)) for x in score_list_h1]
    score_list_h3_sorted = [(x[0], sorted(x[1], reverse=True)) for x in score_list_h3]
    score_list_n1_sorted = [(x[0], sorted(x[1], reverse=True)) for x in score_list_n1]
    score_list_n2_sorted = [(x[0], sorted(x[1], reverse=True)) for x in score_list_n2]


    fig = plt.figure(figsize=(16, 8))

    ax1 = fig.add_subplot(2, 4, 1)
    for k in range(0, 8):
        ax1.plot(score_list_h1[k][1], label=score_list_h1[k][0])
    plt.title("The coverage of H1")  
    
    ax2 = fig.add_subplot(2, 4, 2)
    for k in range(0, 8):
        ax2.plot(score_list_h3[k][1], label=score_list_h3[k][0])
    plt.title("The coverage of H3")

    ax3 = fig.add_subplot(2, 4, 3)
    for k in range(0, 8):
        ax3.plot(score_list_n1[k][1], label=score_list_n1[k][0])
    plt.title("The coverage of N1")

    ax4 = fig.add_subplot(2, 4, 4)
    for k in range(0, 8):
        ax4.plot(score_list_n2[k][1], label=score_list_n2[k][0])
    plt.title("The coverage of N2")

    ax5 = fig.add_subplot(2, 4, 5)
    for k in range(0, 8):
        ax5.plot(score_list_h1_sorted[k][1], label=score_list_h1_sorted[k][0])
    
    ax6 = fig.add_subplot(2, 4, 6)
    for k in range(0, 8):
        ax6.plot(score_list_h3_sorted[k][1], label=score_list_h3_sorted[k][0])

    ax7 = fig.add_subplot(2, 4, 7)
    for k in range(0, 8):
        ax7.plot(score_list_n1_sorted[k][1], label=score_list_n1_sorted[k][0])

    ax8 = fig.add_subplot(2, 4, 8)
    for k in range(0, 8):
        ax8.plot(score_list_n2_sorted[k][1], label=score_list_n2_sorted[k][0])

    plt.legend(bbox_to_anchor=(1.3, 0.5), loc='right', fontsize='small')
    fig.savefig('test3.jpg')

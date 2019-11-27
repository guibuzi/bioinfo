import re
import os
import time
from collections import Counter
from multiprocessing import Pool, Queue, Process, Manager

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


def hamming_dis(a, b):
    count = 0
    for i, j in zip(a, b):
        if i == j:
            count += 1
    return len(a) - count


def cal_most_coverage_aa(epitope_fragments, k):
    """
    计算给定表位集中，存在k个突变的情况下的coverage
    paramters： epitope set, k
    return: coverage
    """
    dis_matrix = []
    score_matrix = []
    for i in range(len(epitope_fragments)):
        dis_row = []
        for j in range(len(epitope_fragments)):
            dis = hamming_dis(epitope_fragments[i], epitope_fragments[j])
            dis_row.append(dis)
        dis_matrix.append(dis_row)
        score_matrix.append([dis_row.count(i) for i in range(9)])
    
    best_score = 0
    best_epitope = ''
    for epitope, score in zip(epitope_fragments, score_matrix):
        if sum(score[:k+1]) > best_score:
            best_epitope = epitope
            best_score = sum(score[:k+1])
    return best_epitope, best_score/len(epitope_fragments)


def cal_aa_pct(seq_list, k, q):
    print("Run task %s" %k)
    length = len(seq_list[0])
    pct_list = []
    for i in range(0, length-10):
        epitope_fragments = [seq[i: i+10] for seq in seq_list]
        aa_count = cal_most_coverage_aa(epitope_fragments, k)
        pct_list.append(aa_count[1])
    print("Task %s done." %k)
    q.put((k, pct_list))



if __name__ == '__main__':
    print("Start task ... ")
    start = time.time()
    data = read_fasta("mosaic/us_ha_sampling_align.fasta")
    seq_list = [value for _, value in data.items()]


    # cal_aa_pct(seq_list, 0)
    # plt.plot(pct_k[0][1])
    # plt.savefig('a3.jpg')

    # pct_k = [cal_aa_pct(seq_list, i) for i in range(9)]
    # for i in range(9):   
    #     plt.plot(pct_k[i])
    # plt.savefig('a3.jpg')

    manager = Manager()
    q = manager.Queue()
    p = Pool(9)
    for i in range(9):
        p.apply_async(cal_aa_pct, args=(seq_list, i, q))
    print("Waiting for all subprocesses done ...")
    p.close()
    p.join()

    values = [q.get() for i in range(9)]
    print("All subprocesses done.")

    values.sort(key=lambda x: x[0])
    for k in values:
        plt.plot(k[1], label=k[0])
    plt.legend()
    plt.savefig('a4.jpg')

    end = time.time()
    print('Cost %6.2f seconds'%(end-start))

# -*- coding: utf-8 -*-

import sys, argparse
import re
import random
from collections import Counter
from pprint import pprint
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


def find_crossover_point(seq1, seq2):
    len1 = len(seq1)
    len2 = len(seq2)
    strings_of_seq1 = set([seq1[i:i+8] for i in range(len1-7)])
    strings_of_seq2 = set([seq2[i:i+8] for i in range(len2-7)])
    common_strings = strings_of_seq1.intersection(strings_of_seq2)
    if len(common_strings) >= 2:
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


def generate_population(seq_list, k, n):
    """
    seq_list: input nature sequences
    k: number of populations
    n: population sizes
    """
    populations = []
    for _ in range(k):
        population_i = []
        popu_size = 0
        while popu_size < n:
            seq1 = random.choice(seq_list)
            seq2 = random.choice(seq_list)
            locations = find_crossover_point(seq1, seq2)
            if locations:
                child1, child2 = recombine(seq1, seq2, locations)
                population_i.append(child1)
                population_i.append(child2)
                popu_size += 1
            else:
                continue
        populations.append(population_i)
    return populations


def cal_coverage(to_cal):
    """
    to_cal: the mosaic to evaluate (list of sequences)
    test_set: nature epitope (list of epitope list)
    return: the coverage of input mosaic
    """
    epitope_mosaic = set([seq[i:i+epitope_length] for seq in to_cal for i in range(len(seq) + 1 - epitope_length)])
    covers = [epitope_count.get(epi, 0) for epi in epitope_mosaic]
    return sum(covers) / num_epitope_nature


def cal_fitness(to_cal, mosaic, k):
    """
    to_cal: the sequence to calculate (string)
    mosaic: current mosaic (list of sequences)
    test_set: nature epitope (list of epitope list)
    k: the index of the population generate this sequence
    return: fitness
    """
    tmp = mosaic[:]
    tmp.pop(k)
    tmp.insert(k,to_cal)
    fitness = cal_coverage(tmp)
    return fitness


def generate_one_parent(popu_i, k):
    seq1 = random.choice(popu_i)
    seq2 = random.choice(popu_i)
    fitness1 = cal_fitness(seq1, mosaic, k)
    fitness2 = cal_fitness(seq2, mosaic, k)
    seq = seq1 if fitness1 >= fitness2 else seq2
    return seq


def if_raw_epitope(children):
    want = []
    for child in children:
        child_epitope = set([child[i:i+9] for i in range(len(child)+1-epitope_length)])
        if child_epitope.issubset(no_raw_epitope):
            want.append(child)
    return want


def generate_child(popu_i, i, inter_crossover_prob=0.5):
    parent_1 = generate_one_parent(popu_i, i)
    if random.random() < inter_crossover_prob:
        parent_2 = generate_one_parent(popu_i, i)  # to do: random selection from nature sequences
    else:
        parent_2 = random.choice((popu_i))
    #print("parent1: ", parent_1)
    #print("parent2: ", parent_2)
    crossover_locations = find_crossover_point(parent_1, parent_2)
    #print("\ncrossover locations:", crossover_locations, "\n")
    if crossover_locations:
        children = recombine(parent_1, parent_2, crossover_locations)
        return children
    else:
        return generate_child(popu_i, i)



if __name__ == "__main__":
    # paramters = sys.argv
    # in_file = paramters[1]
    # log_file = paramters[2]
    # coverage_his_pic = paramters[3]
    # paramters = ["","all_h3.fasta", "logtest.txt", "test"]

    parser = argparse.ArgumentParser(description='Paramters Description')
    parser.add_argument('--in_file', '-i', help='输入序列，FASTA格式 必要参数', required=True)
    parser.add_argument('--log', '-l', help='日志 必要参数，但是有默认值', required=True)
    parser.add_argument('--history', '-his', help='迭代图，必要参数', required=True)
    parser.add_argument('--inter_crossover_prob', '-icp', help='重组概率 属性，必要参数', default=0.5)
    parser.add_argument('--epitope_length', '-e', help='表位长度 属性，必要参数，但是有默认值', default=9)
    parser.add_argument('--raw_threshold', '-r', help='raw_threshold 属性，必要参数，但是有默认值', default=3)
    parser.add_argument('--population_number', '-n', help='population_number 属性，必要参数，但是有默认值', default=4)
    parser.add_argument('--population_size', '-s', help='population_size 属性，必要参数，但是有默认值', default=500)
    args = parser.parse_args()

    in_file = args.in_file
    log_file = args.log
    coverage_his_pic = args.history

    epitope_length = int(args.epitope_length)
    raw_threshold = int(args.raw_threshold)
    inter_crossover_prob = int(args.inter_crossover_prob)
    population_number = int(args.population_number)
    population_size = int(args.population_size)

    # 读入序列 并产生天然表位
    data = read_fasta("/home/zeng/python_work/bioinfo/mosaic/{}".format(in_file))
    seq_list = list(set(data.values()))
    num_nature_seqs = len(seq_list)
    epitope_nature = [seq[i:i+epitope_length] for seq in seq_list for i in range(len(seq)-epitope_length)]
    num_epitope_nature = len(epitope_nature)
    epitope_count = Counter(epitope_nature)
    no_raw_epitope = {k for k, v in epitope_count.items() if v > raw_threshold}   # 超参数 raw number

    # 产生种群 初始 mosaic
    popus = generate_population(seq_list, population_number, population_size)
    init_mosaic = [random.choice(popus[i]) for i in range(population_number)]
    init_coverage = cal_coverage(init_mosaic)


    print("There are {} no redundant sequences in the input file.".format(num_nature_seqs))
    print("Generate %s populations." % population_number)
    print("The size of populations is: %s." % population_size)
    print("The epitope length is: %s" % epitope_length)
    print("Raw epitope threshold is: %s\n" % raw_threshold)
    print("Initial mosaic: ")
    for mm in init_mosaic:
        print(mm)
    print("Initial coverage: {}\n\n".format(init_coverage))

    initlog = open("/home/zeng/Desktop/{}".format(log_file), "w")
    initlog.write("There are {} no redundant sequences in the input file.\n".format(num_nature_seqs))
    initlog.write("Generate %s populations.\n" % population_number)
    initlog.write("The size of populations is: %s.\n" % population_size)
    initlog.write("The epitope length is: %s\n" % epitope_length)
    initlog.write("Raw epitope threshold is: %s\n\n" % raw_threshold)
    initlog.write("Initial mosaic:\n")
    for mm in init_mosaic:
        initlog.write(mm + "\n")
    initlog.write("\nInitial coverage: {}\n\n".format(init_coverage))
    initlog.close()


    # 进入主程序 初始化主程序参数
    mosaic = init_mosaic
    coverage = init_coverage
    old_restart_coverage = init_coverage
    iters = 0
    count_0 = 0
    count_restart_without_improve = 0
    coverage_list = [coverage]
    while True:  # 无限迭代
        log = open("/home/zeng/Desktop/{}".format(log_file), "a")
        print("Iteration {}...\n".format(iters))
        old_coverage = coverage

        # 每一轮迭代中 依次迭代种群
        for i, popu_i in enumerate(popus):
            count = 0
            print("Start iterate {} population.".format(i))
            # 每个种群内部 产生 10 个子代后进入下一种群
            while count < 10:
                # if count == 0:
                #     print("Start new iteration [{}] in population {}: ...".format(count, i))
                children = generate_child(popu_i, i)
                no_raw_children = if_raw_epitope(children)
                count += 1
                if no_raw_children:
                    child = random.choice(no_raw_children)
                    child_score = cal_fitness(child, mosaic, i)
                    to_compare_index = random.sample(list(range(len(popu_i))), 4)
                    to_compare = [popu_i[i] for i in to_compare_index]
                    to_compare_score = [cal_fitness(seq, mosaic, i) for seq in to_compare]
                    #print("child: ", child)
                    #print("child fitness: ", child_score)
                    #print("\nrandom selected 4 sequences's fitness:")
                    #print(to_compare_score)

                    # 如果 child 由于群体，更新群体
                    for index, to_compare in zip(to_compare_index, to_compare_score):
                        if child_score > to_compare:
                            popu_i.pop(index)
                            popu_i.insert(index, child)
                            #print("Update {} to population {}".format(child, i))

                    # 如果 child 优于 mosaic，更新 mosaic
                    if child_score > coverage:
                        mosaic.pop(i)
                        mosaic.insert(i, child)
                        coverage = child_score
                        print("\nUpdate {} to mosaic {}.".format(child, i))
                        #pprint("Current mosaic: ")
                        #pprint(mosaic)
                        print("Current coverage: {}\n".format(coverage))
                else:
                    #print("Reject due to raw epitope")
                    pass

        if iters % 100 == 0:
            log.write("After %s iterations, the current mosaic is:\n" % iters)
            for mm in mosaic:
                log.write(mm + "\n")
            log.write("The current coverage is: %s\n\n" % coverage)

        iters += 1
        coverage_list.append(coverage)
        delta = coverage - old_coverage
        print("\nAfter this iteration current coverage is: {}.".format(coverage))
        print("The fitness have been improved {} during this iteration.\n".format(delta))

        # 判断是否需要产生新的群体
        if delta == 0:
            count_0 += 1
            if count_0 > 50:
                count_0 = 0
                popus = generate_population(seq_list, population_number, population_size)
                delta_between_restart = coverage - old_restart_coverage
                old_restart_coverage = coverage
                if delta_between_restart == 0:
                    count_restart_without_improve += 1
                print("\nRestart populations...\n\n")

        if iters == 1000000 or count_restart_without_improve >= 500:
            print("\nTotal cost {} itearations.".format(iters))
            print("Finanl mosaic:")
            for mm in mosaic:
                print(mm)
            break

        log.close()

    plt.title("Coverage (%) v.s. iteration")
    plt.plot(coverage_list)
    plt.savefig("/home/zeng/Desktop/coverage_history_{}.jpg".format(coverage_his_pic))

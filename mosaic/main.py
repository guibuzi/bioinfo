import sys
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
    strings_of_seq1 = set([seq1[i:i+8] for i in range(len1-8)])
    strings_of_seq2 = set([seq2[i:i+8] for i in range(len2-8)])
    common_strings = strings_of_seq1.intersection(strings_of_seq2)
    if common_strings:
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


def cal_coverage(to_cal, test_set):
    """
    to_cal: the mosaic to evaluate (list of sequences)
    test_set: nature epitope (list of epitope list)
    return: the coverage of input mosaic
    """
    epitope_mosaic = set([seq[i:i+epitope_length] for seq in to_cal for i in range(len(seq)-epitope_length)])
    covers = [epi for epi in epitope_nature if epi in epitope_mosaic]
    return len(covers) / len(epitope_nature)

def cal_fitness(to_cal, mosaic, test_set, k):
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
    fitness = cal_coverage(tmp, test_set)
    return fitness


def generate_one_parent(popu_i, k):
    seq1 = random.choice(popu_i)
    seq2 = random.choice(popu_i)
    fitness1 = cal_fitness(seq1, mosaic, test_set, k)
    fitness2 = cal_fitness(seq2, mosaic, test_set, k)
    seq = seq1 if fitness1 >= fitness2 else seq2
    return seq


def if_raw_epitope(children):
    want = []
    for child in children:
        child_epitope = set([child[i:i+9] for i in range(len(child)-epitope_length)])
        if child_epitope.issubset(no_raw_epitope):
            want.append(child)
    return want


def generate_child(popu_i, i):
    parent_1 = generate_one_parent(popu_i, i)
    parent_2 = generate_one_parent(popu_i, i)  # to do: random selection from nature sequences
    #print("parent1: ", parent_1)
    #print("parent2: ", parent_2)
    crossover_locations = find_crossover_point(parent_1, parent_2)
    #print("\ncrossover locations:", crossover_locations, "\n")
    children = recombine(parent_1, parent_2, crossover_locations)
    return children



if __name__ == "__main__":
    paramters = sys.argv
#    paramters = ["","all_h3.fasta", "logtest.txt"]
    epitope_length = 9
    raw_threshold = 3
    population_number = 4
    population_size = 500

    # 读入序列 并产生天然表位
    data = read_fasta("/home/zeng/python_work/bioinfo/mosaic/{}".format(paramters[1]))
    seq_list = list(set(data.values()))
    epitope_nature = [seq[i:i+epitope_length] for seq in seq_list for i in range(len(seq)-epitope_length)]
    test_set = epitope_nature
    # test_set = [[seq[i:i+epitope_length] for i in range(len(seq)-epitope_length)] for seq in seq_list]
    # nr_epitope_nature = set(epitope_nature)
    epitope_count = Counter(epitope_nature)
    no_raw_epitope = {k for k, v in epitope_count.items() if v > raw_threshold}   # 超参数 raw number

    # 产生种群 初始 mosaic
    popus = generate_population(seq_list, population_number, population_size)
    init_mosaic = [random.choice(popus[i]) for i in range(population_number)]
    init_coverage = cal_coverage(init_mosaic, test_set)


    print("There are {} no redundant sequences in the input file.".format(len(seq_list)))
    print("Generate %s populations." % population_number)
    print("The size of populations is: %s." % population_size)
    print("The epitope length is: %s" % epitope_length)
    print("Raw epitope threshold is: %s\n" % raw_threshold)
    print("Initial mosaic: ")
    for mm in init_mosaic:
        print(mm)
    print("Initial coverage: {}\n\n".format(init_coverage))

    initlog = open("/home/zeng/Desktop/{}".format(paramters[2]), "w")
    initlog.write("There are {} no redundant sequences in the input file.\n".format(len(seq_list)))
    initlog.write("Generate %s populations.\n" % population_number)
    initlog.write("The size of populations is: %s.\n" % population_size)
    initlog.write("The epitope length is: %s\n" % epitope_length)
    initlog.write("Raw epitope threshold is\n\n: %s" % raw_threshold)
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
        log = open("/home/zeng/Desktop/{}".format(paramters[2]), "a")
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
                    child_score = cal_fitness(child, mosaic, test_set, i)
                    to_compare_index = random.sample(list(range(len(popu_i))), 4)
                    to_compare = [popu_i[i] for i in to_compare_index]
                    to_compare_score = [cal_fitness(seq, mosaic, test_set, i) for seq in to_compare]
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
                        print("Update {} to mosaic {}.".format(child, i))
                        #pprint("Current mosaic: ")
                        #pprint(mosaic)
                        print("Current coverage: {}".format(coverage))
                else:
                    #print("Reject due to raw epitope")
                    pass

        if iters % 100 == 0:
            log.write("After %s iterations, the current mosaic is:\n")
            for mm in mosaic:
                log.write(mm + "\n")
            log.write("The current coverage is: %s\n\n" % coverage)

        delta = coverage - old_coverage
        iters += 1
        coverage_list.append(coverage)
        print("\nAfter this iteration current coverage is: {}.".format(coverage))
        print("The fitness have been improved {} during this iteration.\n".format(delta))

        # 判断是否需要产生新的群体
        if delta == 0:
            count_0 += 1
            if count_0 > 50:
                count_0 = 0
                popus = generate_population(seq_list, population_number, len(seq_list))
                delta_between_restart = coverage - old_restart_coverage
                old_restart_coverage = coverage
                if delta_between_restart == 0:
                    count_restart_without_improve += 1
                print("\nRestart populations...\n\n")

        if iters == 100000 or count_restart_without_improve >= 50:
            print("\n\nTotal cost {} itearations.".format(iters))
            print("Finanl mosaic:")
            for mm in mosaic:
                print(mm)
            break

        log.close()

    plt.title("Coverage (%) v.s. iteration")
    plt.plot(coverage_list)
    plt.savefig("/home/zeng/Desktop/coverage_history.jpg")

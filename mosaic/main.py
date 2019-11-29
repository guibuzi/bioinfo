import sys
import s01_initialize as init
import s03_score_function as cal
import random
from collections import Counter
from pprint import pprint


def generate_one_parebt(popu_i, i):
    seq1 = random.choice(popu_i)
    seq2 = random.choice(popu_i)
    fitness1 = cal.cal_fitness(seq1, mosaic, seq_list, i)
    fitness2 = cal.cal_fitness(seq2, mosaic, seq_list, i)
    seq = seq1 if fitness1 >= fitness2 else seq2
    return seq



if __name__ == "__main__":
    data = init.read_fasta("/home/zeng/python_work/bioinfo/mosaic/us_ha_sampling.fasta")
    seq_list = list(data.values())
    epitope_nature = [seq[i:i+9] for seq in seq_list for i in range(len(seq)-9)]
    epitope_count = Counter(epitope_nature)
    no_raw_epitope = {k for k, v in epitope_count.items() if v > 3}

    k = 4
    popus = init.main(seq_list, k, len(seq_list))
    init_mosaic = [random.choice(popus[i]) for i in range(k)]
    init_coverage = cal.cal_coverage(init_mosaic, seq_list)

    mosaic = init_mosaic
    coverage = init_coverage
    print("Initial mosaic: ")
    pprint(mosaic)
    print("\nInitial coverage: {}\n\n\n".format(coverage))

    initlog = open("/home/zeng/Desktop/log.txt", "w")
    initlog.close()
    log = open("/home/zeng/Desktop/log.txt", "a")
    log.write("Initial mosaic: ")
    for mm in mosaic:
        log.write(mm)
    log.write("\nInitial coverage: {}\n\n".format(coverage))


    delta = coverage
    iters = 0
    coverage_list = [coverage]
    while True:
        log.write("The {} iteration.".format(iters))
        print("The {} iteration.".format(iters))
        old_coverage = coverage

        for i, popu_i in enumerate(popus):
            #print("Start iterate {} population.\n".format(i))
            count = 0
            while count < 11:
                count += 1
                #print("Start new iteration [{}] in population {}: ...".format(count, i))
                parent_1 = generate_one_parebt(popu_i, i)
                parent_2 = generate_one_parebt(popu_i, i)
                #print("parent1: ", parent_1)
                #print("parent2: ", parent_2)

                crossover_locations = init.find_crossover_point(parent_1, parent_2)
                #print("\ncrossover locations:", crossover_locations, "\n")
                children = init.recombine(parent_1, parent_2, crossover_locations)
                want = []
                for child in children:
                    child_epitope = set([child[i:i+9] for i in range(len(child)-9)])
                    
                    if child_epitope.issubset(no_raw_epitope):
                        want.append(child)

                if want:
                    child = random.choice(want)
                    #print("child: ", child)
                    child_score = cal.cal_fitness(child, mosaic, seq_list, i)
                    to_compare_index = random.sample(list(range(len(popu_i))), 4)
                    to_compare = [popu_i[i] for i in to_compare_index]
                    to_compare_score = [cal.cal_fitness(seq, mosaic, seq_list, i) for seq in to_compare]
                    #print("child fitness: ", child_score)
                    #print("\nrandom selected 4 sequences's fitness:")
                    #print(to_compare_score)

                    for index, to_compare in zip(to_compare_index, to_compare_score):
                        if child_score > to_compare:
                            popu_i.pop(index)
                            popu_i.insert(index, child)
                            #print("Update {} to population {}".format(child, i))
                    
                    if child_score > coverage:
                        mosaic.pop(i)
                        mosaic.insert(i, child)
                        coverage = child_score
                        log.write("Update {} to mosaic {}.".format(child, i))
                        print("\nUpdate {} to mosaic {}.".format(child, i))
                        #pprint("Current mosaic: ")
                        #pprint(mosaic)
                        log.write("Current coverage: {}".format(coverage))
                        print("Current coverage: {}".format(coverage))
                else:
                    #print("Rejected because raw epitope!")
                    pass
                #print("\n\n")
            
        delta = coverage - old_coverage
        iters += 1
        coverage_list.append(coverage)
        log.write("The fitness have been improved {} during this iteration.\n".format(delta))
        print("The fitness have been improved {} during this iteration.\n".format(delta))
    print("Total cost {} itearations.".format(iters))
    log.close()
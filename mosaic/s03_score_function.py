def cal_coverage(to_cal, nature_seqs):
    epitope_nature = set([seq[i:i+9] for seq in nature_seqs for i in range(0, len(seq)-9)])
    epitope_mosaic = set([seq[i:i+9] for seq in to_cal for i in range(0, len(seq)-9)])
    return len(epitope_mosaic) / len(epitope_nature)


def cal_fitness(to_cal, mosaic, nature_seqs, k):
    tmp = mosaic[:]
    tmp.pop(k)
    tmp.insert(k,to_cal)
    new_coverage = cal_coverage(tmp, nature_seqs)
    return new_coverage

def read_fasta(_in):
    data = {}
    with open(_in) as f:
        for line in f:
            if line.startswith('>'):
                label = line[1:-1]
                data[label] = ''
            else:
                data[label] += line.strip()
    return data

fas_ = read_fasta("/home/zeng/Desktop/riskModel/SIV-NCBI-H1.align.fasta")


def gap_pect(sequence):
    return list(sequence).count('-') / len(sequence)


f = open("/home/zeng/Desktop/riskModel/SIV-NCBI-H1.align.unify.fasta", "w")
for k, v in fas_.items():
    if gap_pect(v) < 0.2:
        f.write(">%s\n%s\n" % (k, v))
f.close()
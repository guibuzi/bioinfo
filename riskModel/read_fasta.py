import json
import pprint


path = "/home/zeng/Desktop/riskModel/SIV-NCBI.fasta"

def read_fasta(_in):
    data = {}
    with open(_in) as f:
        for line in f:
            if line.startswith('>'):
                acc, *(name_), date, _, _ = line[1:-1].split()
                data[acc] = ''
                name = ' '.join(name_)
                # print(acc, name, date, sep='\t')
            else:
                data[acc] += line.strip()
    return data

fas_ = read_fasta(path)


def read_info(_in):
    with open(_in) as f:
        contents = json.load(f)
    return contents

info_ = read_info("/home/zeng/Desktop/riskModel/SIV-NCBI-info-new.txt")


fas2_ =  {k: fas_[k] for k in list(info_.keys())}
print(len(fas_))


with open("/home/zeng/Desktop/riskModel/SIV-NCBI-H1.fasta", "w") as f:
    for k in list(info_.keys()):
        f.write(">%s|%s|%s|%s\n%s\n" % (k, info_[k]['name'], info_[k]['country'], info_[k]['date'], fas_[k]))
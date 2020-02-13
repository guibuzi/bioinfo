import sys, os
import multiprocessing
from collections import defaultdict

def read_acs(_in):
    with open(_in) as f:
        s = f.readlines()
    return [x.strip() for x in s]

def retrive(_in):
    command = "esummary -db nuccore -id %s | xtract -pattern DocumentSummary -element Caption,TaxId" % _in
    result = os.popen(command)  
    res = result.read().strip()
    if res:
        print(res)
    else:
        print("%s\tNo result" % _in)
    return [x.split('\t') for x in res.split('\n')]

def batch_retrive(_in, _out):
    print("Begin retrive %s" % _in)
    command = "cat %s | epost -db nuccore | esummary -db nuccore | xtract -pattern DocumentSummary -element Caption,TaxId > %s" % (_in, _out)
    os.system(command)
    print("Finish retrive %s" % _in)


if __name__ == "__main__":
    acs = read_acs("/home/zeng/Desktop/cov_complete_ac.txt")
    acs = acs[2889:]
    results = []
    num = len(acs)
    count = num // 100 + 1
    for i in range(count):
        result = retrive(repr(' '.join(acs[i*100:(i+1)*100])))
        results.extend(result)

    results_to_write = ["\t".join(i) for i in results]
    with open("/home/zeng/Desktop/cov_complete_ac.result.txt", "w") as f:
        for a in results:
            f.writelines(results_to_write)
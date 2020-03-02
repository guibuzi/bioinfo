from collections import defaultdict
import csv


def read_log(_in):
    results = []
    with open(_in) as f:
        for line in f:
            r = line[:-1].split("\t")
            results.append(r)
    return results

def generate_fragments(qstarts, stitles):
    file_header  =['start', 'stop', 'length', 'subject']
    # csv_file = open('out.csv', 'w')
    # writer = csv.writer(csv_file)
    # writer.writerow(file_header)
    fragments = defaultdict(list)
    start = qstarts[0]
    title = stitles[0]
    for i in range(1, len(stitles)):
        if stitles[i] != stitles[i-1]:
            end = qstarts[i-1]
            fragments[title].append((start, end))
            print("%s : %s-%s, length: %s" % (title, start, end, end-start))
            # writer.writerow([start, end, end-start, title])
            start = qstarts[i]
            title = stitles[i]
    fragments[title].append((start, qstarts[-1]))
    print("%s : %s-%s, length: %s" % (title, start, qstarts[-1], qstarts[-1]-start))
    # writer.writerow([start, qstarts[-1], qstarts[-1]-start, title])
    # csv_file.close()
    return fragments



def main(_in):
    b = read_log("/home/zeng/python_work/bioinfo/nCov/result/%s.log" % _in)
    _, qstarts, qends, _, _, saccs, stitles, snames, _, _, pidents = zip(*b)
    qstarts = [int(x) for x in qstarts if x != 'None']
    qends = [int(x) for x in qends if x != 'None']

    saccs = [x for x in saccs if x != 'None']
    pidents = [float(x) for x in pidents if x != 'None']
    midpoints = [(s+e)//2 for s, e in zip(qstarts, qends)]
    stitles = [x for x in stitles if x != 'None']
    snames = [x for x in snames if x != 'None']

    with open("/home/zeng/Desktop/source_of_pangolin_gd.txt", "w") as f:
        f.write("EPI_ISL_410721\n")
        f.write('\n'.join(set(saccs)))


    generate_fragments(midpoints, stitles)


main('pangolin-gd-1000')
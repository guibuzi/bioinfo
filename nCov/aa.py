import matplotlib.pyplot as plt
from collections import defaultdict, Counter


def read_log(_in):
    results = []
    with open(_in) as f:
        for line in f:
            r = line[:-1].split("\t")
            results.append(r)
    return results


def plot_fig(qstarts, pidents, stitles, _out):
    results = defaultdict(list)
    for qstart, pident, title in zip(qstarts, pidents, stitles):
        results[title].append((qstart, pident))

    fig, ax = plt.subplots(figsize=(18, 4))
    for k, v in results.items():
        xs, ys = zip(*v)
        ax.scatter(xs, [105]*len(xs), label=k)
        ax.scatter(xs, ys, s=1)
        plt.ylabel('identity %')

    plt.legend(bbox_to_anchor=(1,1))
    fig.savefig("/home/zeng/python_work/bioinfo/%s.jpg" % _out)


def color_generator(labels):
    counter = Counter(labels)
    counter.most_common(5)


def main(_in):
    b = read_log("/home/zeng/python_work/bioinfo/%s.log" % _in)
    #b = [x[-4]='P' for x in b if x[-4] == 'N/A']
    _, qstarts, qends, _, _, _, stitles, snames, _, _, pidents = zip(*b)
    qstarts = [int(x) for x in qstarts if x != 'None']
    pidents = [float(x) for x in pidents if x != 'None']
    stitles = list(map(lambda  x: 'BetaCoV/pangolin' if x == 'BetaCov/pangolin' else x, stitles))
    stitles = [x for x in stitles if x != 'None']
    snames = [x for x in snames if x != 'None']

    crossover_point = []
    for i in range(1, len(stitles)):
        if stitles[i] != stitles[i-1]:
            crossover_point.append([qstarts[i-1], stitles[i-1]])
            print([qstarts[i-1], stitles[i-1]])


    #plot_fig(qstarts, pidents, stitles, _in)


main('nCov/2019-ncov-100')

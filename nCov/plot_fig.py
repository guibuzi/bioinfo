import matplotlib.pyplot as plt
from collections import defaultdict

def read_log(_in):
    results = []
    with open(_in) as f:
        for line in f:
            r = line.split("\t")
            results.append(r)
    return results



def plot_fig(qstarts, pidents, stitles, _out):
    results = defaultdict(list)
    for qstart, pident, title in zip(qstarts, pidents, stitles):
        results[title].append((qstart, pident))

    fig, ax = plt.subplots(figsize=(16, 4))
    ax.plot(qstarts, pidents)
    for k, v in results.items():
        xs, ys = zip(*v)
        ax.scatter(xs, [105]*len(xs), label=k)
    #ax.set_ylim(-1, 1)
    #ax.yaxis.set_visible(False)
    plt.legend(bbox_to_anchor=(1, 1))
    fig.savefig("/home/zeng/python_work/bioinfo/%s.jpg" % _out)



def main(_in):
    b = read_log("/home/zeng/python_work/bioinfo/%s.log" % _in)
    _, qstarts, qends, _, _, _, stitles, snames, _, _, pidents = zip(*b)
    qstarts = [int(x) for x in qstarts if x != 'None']
    pidents = [float(x) for x in pidents if x != 'None\n']
    stitles = [x for x in stitles if x != 'None']
    plot_fig(qstarts, pidents, stitles, _in)



tasks = [3000, 2000, 1500, 1000, 500, 250, 200, 150, 100]
for task in tasks:
    main('nCov/2019-ncov-%s' % task)
    # main('nCov/RaTG13-%s' % task)
    # main('nCov/pangolin-%s' % task)
    # main('nCov/ZXC21-%s' % task)

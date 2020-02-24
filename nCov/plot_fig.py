import matplotlib.pyplot as plt
from collections import defaultdict, Counter
import numpy as np
import matplotlib.colors as mcolors


def read_log(_in):
    results = []
    with open(_in) as f:
        for line in f:
            r = line[:-1].split("\t")
            results.append(r)
    return results


def generate_fragments(qstarts, stitles, pidents):
    x_list = defaultdict(list)
    y_list = defaultdict(list)
    for qstart, pident, stitle in zip(qstarts, pidents, stitles):
        x_list[stitle].append(qstart)
        y_list[stitle].append(pident)

    fragments = defaultdict(list)
    start = qstarts[0]
    title = stitles[0]
    for i in range(1, len(stitles)):
        if stitles[i] != stitles[i-1]:
            end = qstarts[i-1]
            fragments[title].append((start, end))
            print("%s : %s-%s" % (title, start, end))
            start = qstarts[i]
            title = stitles[i]
    print("%s : %s-%s" % (title, start, qstarts[-1]))
    fragments[title].append((start, int(qstarts[-1])))
    return fragments, x_list, y_list


def plot_fig(qstarts, stitles, pidents, colormap, _in):
    fragments, x_list, y_list = generate_fragments(qstarts, stitles, pidents)

    fig, ax = plt.subplots(figsize=(16, 4))
    fig.subplots_adjust(0.025, 0.1, 0.85, 0.9)
    #ax.yaxis.set_visible(False)
    for k, v in fragments.items():
        hl = defaultdict(list)
        for item in v:
            hl['xmin'].append(item[0])
            hl['xmax'].append(item[1])
            hl['y'].append(105)
        plt.hlines(y=hl['y'], xmin=hl['xmin'], xmax=hl['xmax'], \
                   label=k, color=colormap[k], linewidth=8)
        plt.scatter(x=x_list[k], y=y_list[k] , color=colormap[k], s=2)
    ax.legend(bbox_to_anchor=(1.01, 1), loc='upper left', borderaxespad=0.)
    fig.savefig("/home/zeng/python_work/bioinfo/%s.jpg" % _in)


def main(_in):
    b = read_log("/home/zeng/python_work/bioinfo/%s.log" % _in)
    _, qstarts, _, _, _, _, stitles, snames, _, _, pidents = zip(*b)
    qstarts = [int(x) for x in qstarts if x != 'None']
    pidents = [float(x) for x in pidents if x != 'None']
    stitles = [x for x in stitles if x != 'None']
    snames = [x for x in snames if x != 'None']

    colormap = {}
    colors = list(mcolors.TABLEAU_COLORS.values())
    top_9 = [i[0] for i in sorted(list(Counter(stitles).most_common(9)), key=lambda x: x[1], reverse=True)]
    print(top_9)
    stitles2 = list(map(lambda x: 'other' if x not in top_9 else x, stitles))
    colormap = {top_9[i]: colors[i] for i in range(len(top_9))}
    colormap['other'] = colors[9]
    print(colormap)

    plot_fig(qstarts, stitles2, pidents, colormap, _in)


tasks = [3000, 2000, 1500, 1000, 500, 250, 200, 150, 100]
for task in tasks:
    main('nCov/RaTG13-%s' % task)
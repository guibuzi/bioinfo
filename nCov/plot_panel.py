import json
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from collections import Counter, defaultdict


key_map = {'RaTG13': 'RaTG13', 
'pangolin/GD': 'pangolin/GD', 
'bat_SL_CoVZXC21': 'bat_SL', 
'bat_SL_CoVZC45': 'bat_SL', 
'pangolin/GX': 'pangolin/GX', 
'Longquan_140': 'Longquan_140', 
'NA': 'NA', 
'HKU3_1': 'HKU3', 'HKU3_2': 'HKU3', 'HKU3_4': 'HKU3', 'HKU3_6': 'HKU3', 'HKU3_3': 'HKU3', 'HKU3_7': 'HKU3', 'HKU3_8': 'HKU3', 'HKU3_13': 'HKU3', 'HKU3_5': 'HKU3', 'HKU3_9': 'HKU3', 'HKU3_11': 'HKU3', 'HKU3_12': 'HKU3', 'HKU3_10': 'HKU3'}

labels = ['RaTG13', 'pangolin/GD', 'pangolin/GX', 'bat_SL', 'Longquan_140', 'HKU3', 'SARSr','NA']
color_map = {v: list(mcolors.TABLEAU_COLORS.values())[i] for i,v in enumerate(labels)}

dirpath='/home/zeng/Desktop/recombination-analysis-0328'

def read_data(_in):
    # read data
    with open(_in) as f:
        results = json.loads(f.read())
    results = {int(k): set(map(lambda x: key_map.get(x[1], 'SARSr'), v)) for k, v in results.items() for i in v if i[0]!='I'}
    # unify by source virus
    source_by_virus = defaultdict(list)
    for k, v in results.items():
        for v_i in v:
            source_by_virus[v_i].append(k)
    return source_by_virus


def fragment(hits, inter_length=20):
    fragments = []
    start = hits[0]
    for i in range(len(hits)-1):
        if hits[i+1] > hits[i] + inter_length:
            end = hits[i] + inter_length - 1
            fragments.append((start, end, end-start+1))
            start = hits[i+1]
    else:
        end = hits[-1]+inter_length-1
        fragments.append((start, end, end-start+1))
    return fragments


def plot_hlines(ax, task, min_length=3, **kwargs):
    source_by_virus = read_data('%s/%s' % (dirpath, task))
    fragments_by_virus = {k: fragment(v, **kwargs) for k, v in source_by_virus.items()}

    yvalue = 0
    yvalues = []
    yticks = []
    for k in labels:
        v = fragments_by_virus.get(k)
        if v:
            xmins, xmaxs, _ = zip(*v)
            for s, e, l in v:
                if l > min_length:
                    ax.text(s, yvalue, str(s), fontsize=6)
                    ax.text(e, yvalue, str(e), fontsize=6)
            ax.hlines(y=[yvalue]*len(v), xmin=xmins, xmax=xmaxs, colors=color_map[k], lw=4)
            yvalues.append(yvalue)
            yticks.append(k)
            yvalue -= 0.1
    ax.set_yticks(yvalues)
    ax.set_yticklabels(yticks)


def main(**kwargs):
    # tasks = ['sars-cov-2_80_500_cds.json', 'ratg13_80_500_cds.json', 'pangolin-gd_80_500_cds.json', 'pangolin-gx_80_500_cds.json']
    tasks = ['2019-ncov-cds.json']
    
    cell_width = 16
    cell_height = 4
    margin = 2
    topmargin = 1

    n = len(tasks)
    width = cell_width * 1 + 2 * margin
    height = cell_height * n + 2 * topmargin

    fig, axes = plt.subplots(n, 1, figsize=(width, height))
    fig.subplots_adjust(margin/width, margin/height, (width-margin)/width, (height-topmargin)/height)

    plot_hlines(axes, tasks[0])
    # for i, task in enumerate(tasks):
    #    plot_hlines(axes[i], task, **kwargs)
    # plt.show()
    plt.savefig('tt.jpg')


main(inter_length=3, min_length=100)
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


def fragment(hits, window_length):
    fragments = []
    start = hits[0]
    for i in range(len(hits)-1):
        if hits[i+1] > hits[i] + window_length:
            end = hits[i] + window_length
            fragments.append((start, end))
            start = hits[i+1]
    else:
        fragments.append((start, hits[-1]+window_length))
    return fragments


def read_data(_in, window_length):
    # read data
    with open("/home/zeng/Desktop/recombination-analysis-0312/%s" % _in) as f:
        results = json.loads(f.read())
    results = {k: set(map(lambda x: key_map.get(x, 'SARSr'), v)) for k, v in results.items()}
    # unify by source virus
    source_by_virus = defaultdict(list)
    for k, v in results.items():
        for v_i in v:
            source_by_virus[v_i].append(int(k))
    # unify by fragments
    fragments_by_virus = {k: fragment(v, window_length) for k, v in source_by_virus.items()}
    return source_by_virus, fragments_by_virus


def plot_scatter(ax, window_length, query, source_by_virus):
    # begin plot
    yvalues = []
    yticks = []
    yvalue = len(source_by_virus)
    for label in labels:
        fragments = source_by_virus.get(label)
        if fragments:
            ax.scatter(x=fragments, y=[yvalue]*len(fragments), color=color_map[label], s=3)
            # update y axis
            yvalues.append(yvalue)
            yticks.append(label)
            yvalue -= 1
        ax.set_yticks(yvalues)
        ax.set_yticklabels(yticks)
        ax.set_title("Query: %s, Window size: %s" % (query, window_length))


def plot_hlines(ax, window_length, query, fragments_by_virus):
    # begin plot
    yvalues = []
    yticks = []
    yvalue = len(fragments_by_virus)
    for label in labels:
        fragments = fragments_by_virus.get(label)
        if fragments:
            xmins, xmaxs = zip(*fragments)
            ax.hlines(y=[yvalue]*len(fragments), xmin=xmins, xmax=xmaxs, color=color_map[label], lw=4)
            # update y axis
            yvalues.append(yvalue)
            yticks.append(label)
            yvalue -= 1
        ax.set_yticks(yvalues)
        ax.set_yticklabels(yticks)
        ax.set_title("Query: %s, Window size: %s" % (query, window_length))


def plot_panel(tasks, plotor):
    '''
    in: tasks is 3-tuple (window_length, file_name, fragments_by_virus)
    '''
    cell_width = 16
    cell_height = 4
    margin = 2
    topmargin = 1

    n = len(tasks)
    width = cell_width * 1 + 2 * margin
    height = cell_height * n + 2 * topmargin

    fig, axes = plt.subplots(n, 1, figsize=(width, height))
    fig.subplots_adjust(margin/width, margin/height, (width-margin)/width, (height-topmargin)/height)

    for i, task in enumerate(tasks):
       plotor(axes[i], task[0], task[2], task[3])
    fig.savefig('/home/zeng/Desktop/%s' % task[1])



if __name__ == '__main__':
    file_names = ['sars-cov-2_80_x_cds.json', 'ratg13_80_x_cds.json', 'pangolin-gd_80_x_cds.json', 'pangolin-gx_80_x_cds.json']
    window_length = 500
    querys = ['SARS-CoV-2', 'RaTG13', 'Pangolin/GD', 'Pangolin/GX']
    tasks = []
    for file_name, query in zip(file_names, querys):
        new_file = file_name.replace('_x_', '_%s_' % window_length)
        source_by_virus, fragments_by_virus = read_data(new_file, window_length)
        tasks.append((window_length, 'window_length_5002222.jpg', query, source_by_virus))
    plot_panel(tasks, plot_scatter)

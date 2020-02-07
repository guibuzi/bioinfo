import matplotlib.pyplot as plt
import numpy as np
import random
import re, os, sys
from collections import Counter, defaultdict


def read_log(_in):
    sacs = []
    snames = []
    pidents = []
    qstarts = []
    qends = []
    with open(_in) as f:
        for line in f:
            row = line.split('\t')
            attrs = row[1].split("|")
            sacs.append(attrs[0])
            snames.append(attrs[1])
            pidents.append(float(row[2]))
            qstarts.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[0]))
            qends.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[1]))
    return sacs, snames, pidents, qstarts, qends 


def judge_in_window(window, chip):
    return False if chip[1]<=window[0] or chip[0]>=window[1] else True

def in_window(window, X, y):
    results = []
    for row in zip(X, y):
        chip = row[0][0:2]
        if judge_in_window(window, chip):
            results.append(list(row))
    return results

def get_winners(X, y, len_of_seq=3822, window_size=25):
    winners = []
    for i in range(len_of_seq-window_size):
        window = (i, i+window_size)
        chips_in = in_window(window, X,  y)
        if chips_in:
            contr = Counter(list(zip(*chips_in))[1])
            winner = contr.most_common(1)[0][0]
            prob = contr.most_common(1)[0][1] / sum(contr.values())
            winners.append((winner, prob))
        else:
            continue
    return winners



def main(_in, fragment, len_of_seq=3822, window_size=25, offset=0):
    _, snames, pidents, qstarts, qends = read_log(_in)

    X = np.array(list(zip(qstarts, qends, pidents)))
    y = np.array(list(map(lambda x: x.upper() , snames)))

    winners = get_winners(X, y, len_of_seq=len_of_seq, window_size=window_size)
    wins = list(zip(*winners))[0]

    index_of_win_i = defaultdict(list)
    for i, win in enumerate(wins):
        index_of_win_i[win].append(i)

    fig, ax = plt.subplots(figsize=(16, 4))
    for k, v in index_of_win_i.items():
        ax.plot(np.array(v)+offset, np.zeros(len(v)), label=k)

    ax.set_xlim(offset, len(wins)+offset)
    ax.set_ylim(-1, 1)
    ax.yaxis.set_visible(False)
    ax.set_title('Recombination Analysis of %s Gene' % fragment, loc='left', fontsize=16)
    plt.legend(bbox_to_anchor=(1, 1))

    def draw_arrow(x, y):
        ax.annotate('%s' % str(x + offset), xy=(x+offset, 0), xytext=(x+offset, y), fontsize=8, 
                arrowprops=dict(arrowstyle='->', connectionstyle="arc3"))

    count = 0
    for k, v in index_of_win_i.items():
        s, e = v[0], v[-1]
        if count == 0:
            count += 1
            continue
        draw_arrow(s, 0.8)
        draw_arrow(e, -0.8)
        print(k, s, e)

    fig.savefig('/home/zeng/python_work/bioinfo/nCov/task/consesus-of-%s.jpg' % fragment)


if __name__ == '__main__':
    # main('/home/zeng/python_work/bioinfo/nCov/task/log-s-4000-3-100', 's')
    # main('/home/zeng/python_work/bioinfo/nCov/task/log-e-4000-3-100', 'e')
    # main('/home/zeng/python_work/bioinfo/nCov/task/log-m-4000-3-100', 'm')
    # main('/home/zeng/python_work/bioinfo/nCov/task/log-n-4000-3-100', 'n')
    main('/home/zeng/python_work/bioinfo/nCov/task/log-orf1ab-4000-3-100', 'orf1ab')
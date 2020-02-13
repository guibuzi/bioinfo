import matplotlib.pyplot as plt
import numpy as np
import random
import re, os, sys
from collections import Counter, defaultdict


def read_log(_in):
    snames = []
    pidents = []
    qstarts = []
    qends = []
    with open(_in) as f:
        for line in f:
            row = line.split('\t')
            snames.append(row[2])
            pidents.append(float(row[3]))
            qstarts.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[0]))
            qends.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[1]))
    return snames, pidents, qstarts, qends 


def judge_in_window(window, chip): # window (start, end)  chip (start, end)
    return False if chip[1]<=window[0] or chip[0]>=window[1] else True

def in_window(window, X, y): # X (start, end, pident) y label
    results = []
    for x_i, y_i in zip(X, y):
        chip = x_i[0], x_i[1]
        if judge_in_window(window, chip):
            results.append([x_i, y_i])
    return results

def get_winner(_in, metric='identity'):
    if metric == 'identity':
        winner = sorted(_in, key=lambda x: x[0][2], reverse=True)[0][1]
    elif metric == 'count':
        winner = Counter(list(zip(*_in))[1]).most_common(1)[0][0]
    return winner


def get_winners(X, y, len_of_seq=30000, window_size=25):
    winners = []
    for i in range(len_of_seq-window_size):
        window = (i, i+window_size)
        chips_in = in_window(window, X,  y)
        if chips_in:
            winner = get_winner(chips_in, metric=metrix)
            winners.append((i, winner))
        else:
            continue
    return winners



def main():
    snames, pidents, qstarts, qends = read_log(_path)

    X = np.array(list(zip(qstarts, qends, pidents)))
    y = np.array(list(map(lambda x: x.upper() , snames)))

    winners = get_winners(X, y, len_of_seq=sequence_length, window_size=window_size)

    index_of_win_i = defaultdict(list)
    for i, winner in winners:
        index_of_win_i[winner].append(i)

    fig, ax = plt.subplots(figsize=(16, 4))
    for k, v in index_of_win_i.items():
        ax.scatter(np.array(v)+offset, np.zeros(len(v)), label=k)

    ax.set_xlim(offset, winners[-1][0]+offset)
    ax.set_ylim(-1, 1)
    ax.yaxis.set_visible(False)
    ax.set_title('Recombination Analysis of %s Gene' % _fragment, loc='left', fontsize=16)
    plt.legend(bbox_to_anchor=(1, 1))

    def draw_arrow(x, y):
        ax.annotate('%s' % str(x + offset), xy=(x+offset, 0), xytext=(x+offset, y), fontsize=8, 
                arrowprops=dict(arrowstyle='->', connectionstyle="arc3"))

    for k, v in index_of_win_i.items():
        for i in range(1, len(v)):
            if v[i] - v[i-1] > 1:
                draw_arrow(v[i-1], 0.8)
                print(k, v[i-1])

    fig.savefig('consesus-of-%s.jpg' % _fragment)


if __name__ == '__main__':
#    _, _path, _fragment, window_size = sys.argv
    _path, _fragment, window_size = "nCov/log-ratg13-all-3000", 'all', '50'
    window_size = int(window_size)
    sequence_length = 30000
    offset=0
    metrix='identity'

    main()
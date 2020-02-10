# encoding: utf-8

import sys, os
import random
import re
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

class RandomCut():
    def __init__(self, sequence='', length=0, mode='s', seg_length='3-100'):
        self.sequence = sequence
        self.length = length
        self.mode = mode
        self.range = self.mode_select(mode)
        self.seg_length = self.extract_seg(seg_length)
    
    def __len__(self):
        return self.length
    
    def extract_seg(self, seg_length):
        return tuple(map(int, seg_length.split('-')))
    
    def mode_select(self, mode):
        '''mode in ['s', 'e', 'm', 'n'] or range xx-xx'''
        if mode == 's':
            return (21563, 25385)
        elif mode == 'e':
            return (26245, 26472)
        elif mode == 'm':
            return (26523, 27191)
        elif mode == 'n':
            return (28274, 29533)
        elif mode == 'orf1ab':
            return (266, 21555)
        else:
            return tuple(map(int, mode.split('-')))

    def read_sequence(self, _in):
        seq = ''
        with open(_in) as f:
            for line in f:
                if not line.startswith(">"):
                    seq += line.replace("\n", "")        
        self.sequence = seq[self.range[0]: self.range[1]]
        self.length = len(self.sequence)

    def random_cut(self):
        len_l, len_u = self.seg_length
        len_of_seg = random.randrange(len_l, len_u+1)
        start = random.randrange(0, self.length - len_of_seg)
        end = start + len_of_seg
        return (self.sequence[start: end], (start, end-1))

    def random_cut_write(self):
        segment, _range = self.random_cut()
        with open('segment-%s.fasta' % self.mode, 'w') as f:
            query = ">segment[%s..%s]\n%s\n" % (_range[0], _range[1], segment)
            f.write(query)


class Blastor():
    def __init__(self, _query, _blastdb):
        self.query = _query
        self.blastdb = _blastdb

    def blastn(self, quiet=False):
        command = "blastn -query %s -db %s -task blastn -outfmt '6 qseqid sseqid pident sstart send length' -num_alignments 10" % (self.query, self.blastdb)
        result = os.popen(command)  
        res = result.read()
        if not quiet:
            if res:
                print(res)
            else:
                print("No hits.")
        return res
    
    def extract_result(self, quiet=False):
        res_blast = self.blastn(quiet=quiet)
        if res_blast:
            return [i for i in res_blast.split('\n')[:1] if i]


class Leaderboard():
    def __init__(self):
        pass


def read_log(_in):
    sacs = []
    snames = []
    pidents = []
    pstarts = []
    pends = []
    with open(_in) as f:
        for line in f:
            if line.startswith('segment'):
                row = line.split('\t')
                attrs = row[1].split("|")
                sacs.append(attrs[0])
                snames.append(attrs[1])
                pidents.append(float(row[2]))
                pstarts.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[0]))
                pends.append(int(re.findall(r'\d*\.\.\d*', row[0])[0].split('..')[1]))
    return sacs, snames, pidents, pstarts, pends 


def judge_exist(_in):
    status = os.path.exists(_in)
    if status:
        os.system('rm -rf %s' % _in) 


def randomcolor(num_of_color):
    color_B = sorted(random.sample(range(0, 256), num_of_color))
    color_hex = list(map(lambda x: "#%s00FA" % str(hex(x))[-2:].replace('x', '0').upper(), color_B))
    return color_hex


def plot_colortable(colors, title, emptycols=0):
    cell_width = 212
    cell_height = 22
    swatch_width = 48
    margin = 12
    topmargin = 40

    names = list(colors)
    n = len(names)
    ncols = 4 - emptycols
    nrows = n // ncols + int(n % ncols > 0)

    width = cell_width * 4 + 2 * margin
    height = cell_height * nrows + margin + topmargin
    dpi = 72

    fig, ax = plt.subplots(figsize=(width / dpi, height / dpi), dpi=dpi)
    fig.subplots_adjust(margin/width, margin/height,
                        (width-margin)/width, (height-topmargin)/height)
    ax.set_xlim(0, cell_width * 4)
    ax.set_ylim(cell_height * (nrows-0.5), -cell_height/2.)
    ax.yaxis.set_visible(False)
    ax.xaxis.set_visible(False)
    ax.set_axis_off()
    ax.set_title(title, fontsize=24, loc="left", pad=10)

    for i, name in enumerate(names):
        row = i % nrows
        col = i // nrows
        y = row * cell_height
        swatch_start_x = cell_width * col
        swatch_end_x = cell_width * col + swatch_width
        text_pos_x = cell_width * col + swatch_width + 7
        ax.text(text_pos_x, y, name, fontsize=14,
                horizontalalignment='left',
                verticalalignment='center')
        ax.hlines(y, swatch_start_x, swatch_end_x, color=colors[name], linewidth=18)
    fig.savefig('color-table-of-%s.jpg' % task)
    return fig


# def get_argument():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('')


def path_proc():
    current_path = os.getcwd()
    judge_exist(os.path.join(current_path, 'segment-%s.fasta' % mode))
    # os.chdir(current_path + "/tmp/")
    return current_path


def main():
    current_path = path_proc()

    nCOV = RandomCut(mode=mode, seg_length=seg_length)
    nCOV.read_sequence('%s' % in_file)
    blastor = Blastor('segment-%s.fasta' % mode, '~/blastdb/coronavirus/all-cov')

    def get_sseq_id():
        nCOV.random_cut_write()
        ac_name_pident_s_e  = blastor.extract_result(quiet=True)
        if ac_name_pident_s_e:
            return ac_name_pident_s_e
        else:
            return get_sseq_id()

    results = []
    judge_exist(os.path.join(current_path, 'log-%s' % task))
    with open(os.path.join(current_path, 'log-%s' % task), 'a') as f:
        for _ in range(iter_nums):
            a = get_sseq_id()
            results.extend(a)
            for i in a:
                print(i)
                f.write(i + "\n")
    
    # qseqid sseqid pident sstart send length = list(zip(*results))

    sacs, snames, pidents, pstarts, pends = read_log(os.path.join(current_path, 'log-%s' % task))

    color_sp = randomcolor(len(set(snames)))
    cmap = {n: color_sp[i] for i, n in enumerate(set(snames))}
    cmap['RaTG13'] = '#FF3030'

    plt.figure(figsize=(16, 9))
    for _, sname, pident, pstart, pend in zip(sacs, snames, pidents, pstarts, pends):
        plt.hlines(pident, pstart, pend, label=sname, color=cmap[sname])
    plt.ylabel('% Identity')
    plt.title("%s gene" % mode.upper())
    plt.savefig(os.path.join(current_path, "result-of-%s.jpg" % task))

    plot_colortable(cmap, "Color Table of %s" % task, emptycols=0)

if __name__ == '__main__':
    params = sys.argv[1:]
    task = params[0]
    in_file = params[1]
    mode = params[2]
    iter_nums = int(params[3])
    seg_length = params[4]

    # tmp ='tmp1'
    # log = 'log1'
    # iter_nums = 1000
    main()

    os.system('rm -rf segment-*')
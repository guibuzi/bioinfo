# -*- coding: utf-8 -*-

import os, sys
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import numpy as np
from multiprocessing import Pool, Process


class Blastor():
    def __init__(self, _query,  _db, _mask):
        self.query = _query
        self.db = _db
        self.mask = _mask
        self.query_seq = self._read_query(_query)
        self.query_len = len(self.query_seq)

    def _read_query(self, _in):
        with open(_in)as f:
            seq = ''
            for line in f:
                if not line.startswith(">"):
                    seq += line.strip()
        return seq

    def _blastn(self, query_loc):
        command = "blastn -query %s -query_loc %s -db %s \
                   -num_alignments 10 -task blastn -max_hsps 1\
                   -outfmt '6 qseqid qstart qend sseqid staxid sacc stitle scomname sstart send pident' \
                   -negative_seqidlist %s"  % (self.query, query_loc, self.db, self.mask)
        result = os.popen(command)  
        res = result.read()
        return res

    def one_blast(self, query_loc='all'):
        if query_loc == 'all':
            result = self._blastn('1-%s' % self.query_len)
        else:
            result = self._blastn(query_loc)
        print('The best alignment for %s is %s' % (query_loc, result.split("\n")[0].split("\t")))
        return result.split("\n")[0].split("\t")

    def window_blast(self, log_file, window_size=1000):
        fw = open("%s" % log_file, "w")
        results = []
        for i in range(1, self.query_len - window_size):
            query_loc = "%s-%s" % (i, i+window_size)
            result = self.one_blast(query_loc)
            if result == ['']:
                results.append([None]*11)
                fw.write("\t".join(['None']*11) + "\n")
            else:
                results.append(result)
                fw.write("\t".join(result) + "\n")
        fw.close()
        return results



def main(_in_file, _mask, window_size, outplot=False):
    test = Blastor('/home/zeng/Desktop/%s' % _in_file, 'all-cov/all-cov', '/home/zeng/Desktop/%s' % _mask)
    b = test.window_blast("%s-%s.log" % (_in_file.split('.')[0], window_size), window_size=int(window_size))
    if outplot:
        _, qstarts, qends, _, _, _, stitles, snames, _, _, pidents = zip(*b)
        results = defaultdict(list)
        for qstart, pident, title in zip(qstarts, pidents, stitles):
            if qstart:
                results[title].append((qstart, pident))

        fig, ax = plt.subplots(figsize=(16, 4))
        for k, v in results.items():
            xs, _ = zip(*v)
            ax.scatter(xs, [0]*(xs), label=k)
    #    ax.set_ylim(-1, 1)
    #    ax.yaxis.set_visible(False)
        plt.legend(bbox_to_anchor=(1, 1))
        fig.savefig("%s-%s.jpg" % (_in_file.split('.')[0], window_size))


if __name__ == "__main__":
    _, _in_file, _mask = sys.argv
    # _in_file = 'RaTG13.fasta'
    # _mask = '2019-ncov-ratg.idlist'
    # window_size = 500

    # qaccs, qstarts, qends, sseqids, staxids, saccs, stitles, snames, sstarts, sends, pidents

    process_list = []
    tasks = [3000, 2000, 1500, 1000, 500, 250, 200, 150, 100]
    for task in tasks:
        p = Process(target=main, args=(_in_file, _mask, task, False, ))
        p.start()
        process_list.append(p)
    
    for p in process_list:
        p.join()

#    main(_in_file, _mask, task, 500)
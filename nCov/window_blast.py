# -*- coding: utf-8 -*-

import os, sys
import matplotlib.pyplot as plt
from collections import Counter, defaultdict
import numpy as np


class Blastor():
    def __init__(self, _query,  _db, _mask):
        self.query = _query
        self.db = _db
        self.mask = _mask
        self.query_seq = self._read_query(_query)
        self.query_len = len(self.query_seq)

    def _blastn(self, query_loc):
        command = "blastn -query %s -query_loc %s -db %s -num_alignments 5 \
                   -outfmt '6 qseqid qstart qend sseqid staxid sacc stitle scomname sstart send pident' \
                   -negative_seqidlist %s" % (self.query, query_loc, self.db, self.mask)
        result = os.popen(command)  
        res = result.read()
        return res
    
    def _read_query(self, _in):
        with open(_in)as f:
            seq = ''
            for line in f:
                if not line.startswith(">"):
                    seq += line.strip()
        return seq
    
    def one_blast(self, query_loc='all'):
        if query_loc == 'all':
            result = self._blastn('1-%s' % self.query_len)
        else:
            result = self._blastn(query_loc)
        return result.split("\n")[0].split("\t")

    def window_blast(self, window_size=1000):
        results = []
        for i in range(1, self.query_len - window_size):
            query_loc = "%s-%s" % (i, i+window_size)
            result = self.one_blast(query_loc)
            if result:
                results.append(result)
            else:
                results.append([None]*11)
        return results


if __name__ == "__main__":
    # _, window_size = sys.argv

    window_size = 100
    test = Blastor('/home/zeng/Desktop/2019-ncov.fasta', 'all-cov/all-cov', '/home/zeng/Desktop/2019-ncov.idlist')

    b = test.window_blast(window_size=int(window_size))
    print(b)
    # _, qstarts, qends, sseqids, staxids, saccs, stitles, snames, sstarts, sends, pidents = zip(*b)

    # results = defaultdict(list)
    # for i, stitle in enumerate(stitles):
    #     results[stitle].append(i)

    # fig, ax = plt.subplots(figsize=(16, 4))
    # for k, v in results.items():
    #     ax.scatter(np.array(v), np.zeros(len(v)), label=k)

    # ax.set_ylim(-1, 1)
    # ax.yaxis.set_visible(False)
    # plt.legend(bbox_to_anchor=(1, 1))
    # plt.savefig("test-%s.jpg" % window_size)

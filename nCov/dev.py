import os, sys 
import numpy as np
import re, copy, time, json
from collections import Counter, defaultdict
from multiprocessing import Pool, Process
from multiprocessing.dummy import Pool as ThreadPool


class Sequence:
    def __init__(self, acc=None, title=None, sequence=None, start=None, end=None):
        self.acc = acc
        self.title = title
        self.sequence = sequence
        self.start = start
        self.end = end
    
    def __len__(self):
        return len(self.sequence)
    
    def load_seqs(self, file):
        self.sequence = ''
        with open(file) as f:
            for line in f:
                if line.startswith(">"):
                    self.acc, *title = line[1:-1].split()
                    self.title = ' '.join(title)
                else:
                    self.sequence += line.strip()
            else:
                self.file = file
                self.start = 1
                self.end = len(self.sequence)
    
    def extract_part(self, start, end):
        start = start if start > 0 else 1
        sequence = self.sequence[start-1: end]
        return Sequence(self.acc, self.title, sequence, start, end)
    
    def fasta(self):
        return ">%s %s..%s..%s\n%s\n" % (self.acc, self.start, self.end, self.title, self.sequence)


class Subject(Sequence):
    def __init__(self, acc, title, sequence, start, end, pident=0, score=0):
        super().__init__(acc, title, sequence, start, end)
        self.pident = pident
        self.score = score
    
    def get_bident(self, to_compare, length=1000, db='all-cov/all-cov'):
        start = self.start-length if self.start-length > 0 else 1
        end = self.end+length
        commond ="blastdbcmd -db %s -entry %s -range %s-%s | \
                  blastn -subject %s -outfmt '6 pident' \
                  -max_hsps 1" % (db, self.acc, start, end, to_compare.file)
        try:
            bident = float(os.popen(commond).read().strip())
        except:
            bident = 0
        return bident
    
    def update(self, seq):
        self.sequence = seq

    def __eq__(self, item):
        return self.acc == item.acc


class Result(list):
    """list of Subject class"""
    def __init__(self, _in=list()):
        for i in _in:
            self.append(i)
        
    def top1(self):
        results = sorted(self, key=lambda x: x.score, reverse=True)
        return Result([subject for subject in results if subject.score == results[0].score])
    
    def get_item(self, item):
        try:
            return [value for value in self if value.acc == item.acc][0]
        except:
            raise KeyError('No such item!')
    
    def __contains__(self, item):
        return item.acc in [i.acc for i in self]

    def intersection(self, item):
        return Result([i for i in self if i in item])


class Query(Sequence):
    def __init__(self, acc, title, sequence, start, end):
        super().__init__(acc, title, sequence, start, end)
        self.results = Result()

    def blastn(self, mask, db='all-cov/all-cov'):
        commond = "echo '%s' | blastn -db %s -max_hsps 1 \
                        -outfmt '6 qacc qstart qend sacc sstart send pident score sseq stitle' \
                        -negative_seqidlist %s -task blastn \
                        2> /dev/null" % (self.sequence, db, mask)
        for x in os.popen(commond).readlines():
            _, _, _, sacc, sstart, send, pident, score, sseq, stitle = x.rstrip().split('\t')
            subject = Subject(sacc, stitle, sseq, int(sstart), int(send), float(pident), int(score))
            self.results.append(subject)

    def get_top1(self):
        return self.results.top1()
    
    def get_item(self, item):
        return self.results.get_item(item)


class Triplet:
    def __init__(self, query, background, subject):
        self.query = query
        self.background = background
        self.subject = subject
        self.n_subject = len(subject)
    
    def fasta(self, write=True):
        query_ = self.query.fasta()
        background_ = self.background.fasta()
        subject_ = ''.join([i.fasta() for i in self.subject])
        if write:
            f = open('%s_%s_%s.fasta' % (self.query.acc, self.query.start, self.query.end), 'w')
            f.write(''.join([query_, background_, subject_]))
            f.close()
        return query_, background_, subject_

    def align(self):
        _in = '%s_%s_%s.fasta' % (self.query.acc, self.query.start, self.query.end)
        self.fasta()
        res = os.popen('mafft --auto --quiet %s' % (_in)).read().split('>')[1:]
        os.remove(_in)
        count = 0
        subjects = Result()
        for seq in res:
            label, *nucl = seq.split('\n')
            acc, *attr = label.split()
            (start, end, title), seq = ' '.join(attr).split('..'), ''.join(nucl)
            sequence = Subject(acc, title, seq, start, end)
            if count < 1:
                query = sequence
            elif count < 2:
                background = sequence
            elif count < 2 + self.n_subject:
                subjects.append(sequence)
            count += 1
        self.query, self.subject, self.background = query, subjects, background

    def _identity(self, x, y):
        count = 0
        for x_i, y_i in zip(x.sequence, y.sequence):
            if x_i == y_i:
                count += 1
        return count / len(x)

    def is_subject(self):
        background_pident = self._identity(self.query, self.background)
        s_pident = [self._identity(self.query, s) for s in self.subject]
        # print('back: ', background_pident, 'subject: ', s_pident)
        for pident in s_pident:
            if background_pident < pident:
                return True

    def is_background(self):
        if self.background in self.subject:
            return True


class Bootstrap:
    def __init__(self, triplet):
        self.triplet = triplet
        self._query = np.array(list(triplet.query.sequence))
        self._background = np.array(list(triplet.background.sequence))
        self._subject = [np.array(list(s.sequence)) for s in triplet.subject]
        self._length = len(triplet.query.sequence)

    def sampling(self, sample_size=1000):
        samples = []
        idx = np.random.randint(0, self._length, size=(sample_size, self._length))
        querys = self._query[idx]
        backgrounds = self._background[idx]
        subjects = zip(*[i[idx] for i in self._subject])
        for query, background, subject in zip(querys, backgrounds, subjects):
            tmp = copy.deepcopy(self.triplet)
            tmp.query.update(''.join(query))
            tmp.background.update(''.join(background))
            for i, s in enumerate(subject):
                tmp.subject[i].update(''.join(s))
            samples.append(tmp)
        return samples
            

def get_background(sequence, mask, start, end):
    sequence_background = sequence.extract_part(start, end)
    query_background = Query(sequence_background.acc, sequence_background.title, sequence_background.sequence, sequence_background.start, sequence_background.end)
    query_background.blastn(mask)
    back = query_background.results.top1()
    return back[0]


def one_window(start, *args, **kwargs):
    end = start + window_size - 1
    b_start = start - b_length
    b_end = end + b_length

    sequence_part = sequence.extract_part(start, end)
    query = Query(sequence_part.acc, sequence_part.title, sequence_part.sequence, sequence_part.start, sequence_part.end)
    query.blastn(mask)
    background = query.get_item(get_background(sequence, mask, b_start, b_end))
    subjects = Result([subject for subject in query.get_top1() if subject.get_bident(sequence) > b_ident_threshod])
    
    triplet = Triplet(query, background, subjects)
    if triplet.is_background():
        winner = [triplet.background.title]
        bootstrap_support = None
    else:
        triplet.align()
        bootstrap = Bootstrap(triplet)
        replicates = bootstrap.sampling(bootstrap_value)
        count= 0
        for replicate in replicates:
            if replicate.is_subject():
                count += 1
        else:
            bootstrap_support = count / bootstrap_value
        if bootstrap_support >= bootstrap_support:
            winner = list(set([s.title for s in triplet.subject]))
        else:
            winner = [triplet.background.title]
    print("Position: %s-%s, Background: %s, Winner: %s, Bootstrap support: %s" % (start, end, triplet.background.title, ' '.join(winner), bootstrap_support))        
    return start, winner

def main(start, stop, step):
    # 创建一个工作者线程池
    pool = Pool(10) 
    results = pool.map(one_window, range(start, stop, step))
    #close the pool and wait for the work to finish
    pool.close() 
    pool.join()
    results = {x: y for x, y in results}
    return results

if __name__ == '__main__':
    os.chdir("/home/zeng/Desktop/recombination-analysis-0328")
    file = "ratg13-cds.fasta"
    mask = 'idlist.2019-ncov-ratg'
    sequence = Sequence()
    sequence.load_seqs(file)

    window_size = 500
    step = 3
    b_length = 1000
    b_ident_threshod = 80
    bootstrap_value = 1000
    bootstrap_support = 0.7
    stop = len(sequence.sequence) - window_size

    #one_window(1)
    results = main(1, stop, step)
    with open('ratg13-%s.json' % window_size, 'w') as f:
        f.write(json.dumps(results, indent='\t', sort_keys=True, separators=(',', ': ')))



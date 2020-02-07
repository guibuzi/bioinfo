import sys, os, re, random
from collections import Counter, defaultdict
import  matplotlib.pyplot as plt


def read_seq(_in):
    labels, seqs = [], []
    with open(_in) as f:
        for line in f:
            label, seq = line.split('\t')
            seqs.append(seq.strip())
            labels.append(label)
    return labels, seqs


path = '/Users/jinfeng/Downloads/test.seq'
size = 9

_, sequences = read_seq(path)

def get_k_mers(_in):
    k_mers = []
    seq_len = len(_in)
    for i in range(len(_in) - size):
        k_mers.append(_in[i: i+size])
    return Counter(k_mers)

def get_all_k_mers(_in):
    cout_all = defaultdict(int)
    for seq in _in:
        cout_i = get_k_mers(seq)
        for k, v in cout_i.items():
            cout_all[k] += v
    return cout_all


count_of_all_kmers = get_all_k_mers(sequences)
n_unique_kmers = len(count_of_all_kmers)
n_kmers = sum(list(count_of_all_kmers.values()))
top_10 = sorted(list(count_of_all_kmers.items()), key=lambda x: x[1], reverse=True)[:10]

n_kmers_no_rare = 0
n_unique_kmers_no_rare = 0
for k, v in count_of_all_kmers.items():
    if v >= 3:
        n_kmers_no_rare += v
        n_unique_kmers_no_rare += 1


print("# of input sequences: ", len(sequences))
print("# of unique kmers: ", n_unique_kmers)
print("# of all kmers: ", n_kmers)
print("# of unique no rare kmers: ", n_unique_kmers_no_rare)
print("# of no rare kmers: ", n_kmers_no_rare)
print("top_10 kmers: ", top_10)


to_evaluate = [
"MGGKWSKSSIVGWPAIRERMRRTEPRTEPAAEGVGAVSQDLARHGAITSSNTAANNPDCAWLEAQEEDE",
"MGGKWSKSSVVGWPAVRERMRRTEPAAEGVGAASQDLDKHGAITSSNTAATNADCAWLEAQEDEE",
"MGGKWSKSSIVGWPAVRERIRRTEPAAEGVGAASRDLEKHGAITSSNTATNNADCAWLQAQEEEE",
"MGSKWSKSSIVGWPAVRERMRRAEPAADGVGAVSRDLERHGAITSSNTAANNADCAWLEAQEEEE"]

# to_evaluate = [
# 'MGGKWSKSSIVGWPAIRERMRRAEPAAEGVGAVSQDLARHGAITSSNTAANNADCAWLQAQEEEE',
# 'MGGKWSKSSIVGWPAIRERMRRTEPAAEGVGAASQDLDKHGAITSSNTATNNADCAWLEAQEEEE',
# 'MGGKWSKSSVVGWPAVRERMRRAEPAADGVGAVSRDLERHGAITSSNTAATNADCAWLEAQEEDE',
# 'MGSKWSKSSIVGWPAVRERIRRTEPAAEGVGAASRDLEKHGAITSSNTAANNPDCAWLEAQEDEE']


kmers_of_query = get_all_k_mers(to_evaluate)
n_unique_kmers_in_query = len(kmers_of_query)
n_kmers_in_query = sum(list(kmers_of_query.values()))

print("# of query sequences: ", len(to_evaluate))
print("# of unique kmers: ", n_unique_kmers_in_query)
print("# of all kmers: ", n_kmers_in_query)

score = 0
score_list = []
for kmer in kmers_of_query:
    score_list.append(count_of_all_kmers.get(kmer, 0))
    score += count_of_all_kmers.get(kmer, 0)

print("lowest score: ", sorted(score_list)[:10])
print("Hightest score: ", sorted(score_list, reverse=True)[:10])
print("score of query: ", score)


print(score / n_kmers)
print(score / n_kmers_no_rare)
print()
print(n_unique_kmers_in_query / n_unique_kmers_no_rare)

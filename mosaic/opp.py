epi_len = 9
raw_threshold = 3
population_number = 4
population_size = 500

def read_fasta(path):
    seq_dict = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                label = line[1: -1]
                seq_dict[label] = ''
            else:
                seq = line.replace("\n", "").strip()
                seq_dict[label] += seq
    return seq_dict

data = read_fasta("/home/zeng/python_work/bioinfo/mosaic/us_ha_sampling.fasta")
seq_list = list(set(data.values()))


class Sequences():
    def __init__(self, seq):
        self.seq = seq
        self.epitope = self.get_epitope()

    def get_epitope(self):
        return [self.seq[i:i+epi_len] for i in range(len(self.seq) - epi_len)]
        

class Mosaic():
    def __init__(self, seqs):
        if not isinstance(seqs, list):
            raise TypeError("Expect list type!")
        if seqs == None:
            raise ValueError("At least one sequence!")
        self.seqs = seqs
        self.size = len(seqs)
    
    def get_coverage(self):
        epitope_mosaic = [seq[i:i+epi_len] for seq in self.seqs
                                           for i in range(len(seq)-epi_len)]
        covers = [epi for epi in epitope_nature if epi in epitope_mosaic]
        return len(covers) / len(epitope_nature)
    
    def update(self, seq, i):
        self.seqs.pop(i)
        self.seqs.insert(i, seq)
    

class Population():
    def __init__(self, seqs):
        self.population = seqs


seq1 = Sequences("MGGKWSKSSIVGWPAIRERMRRTEPRTEPAAEGVGAVSQDLARHGAITSSNTAANNPDCAWLEAQEEDE")
print(seq1.seq)
print(seq1.epitope)

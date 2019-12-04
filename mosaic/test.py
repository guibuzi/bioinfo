from statistics import mean
import random
from collections import Counter

seq_list = []
with open("/home/zeng/Desktop/test.seq") as f:
    for line in f:
        seq = line.split()[1].strip()
        seq_list.append(seq)
epitope_length = 9

mosaic = ["MGGKWSKSSIVGWPAIRERMRRTEPRTEPAAEGVGAVSQDLARHGAITSSNTAANNPDCAWLEAQEEDE",
          "MGGKWSKSSVVGWPAVRERMRRTEPAAEGVGAASQDLDKHGAITSSNTAATNADCAWLEAQEDEE",
          "MGGKWSKSSIVGWPAVRERIRRTEPAAEGVGAASRDLEKHGAITSSNTATNNADCAWLQAQEEEE",
          "MGSKWSKSSIVGWPAVRERMRRAEPAADGVGAVSRDLERHGAITSSNTAANNADCAWLEAQEEEE"]

mosaic2 = ["MKTIIAFSCILCLIFAQKLPGSDNSMATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGRICNSPHQILDGKNCTLIDALLGDPHCDDFQNKEWDLFVERSTAYSNCYPYYVPDYATLRSLVASSGNLEFTQESFNWTGVAQDGSSYACRRGSVNSFFSRLNWLYNLNYKYPEQNVTMPNNDKFDKLYIWGVHHPGTDKDQTNLYVQASGRVIVSTKRSQQAVIPNIGSRPWVRGVSSIISIYWTIVKPGDILLINSTGNLIAPRGYFKIQSGKSSIMRSDAHIDECNSECITPNGSISNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPERQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAINQITGKLNRVIKKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMSKLFERTRRQLRENAEDMGNGCFKIYHKCDNTCIGSIRNGTYDHDIYRNEALNNRFQIKGVQLKSGYKDWILWISFAISCFLLCVVLLGFIMWACQKG",
"MKTIIALSYILCLVFAQKIPGNDNSTATLCLGHHAVPNGTIVKTITSDRIEVTNATELVQNSSIGEICDSPHQILDGGNCTLIDALLGDPQCDGFQNRKWDLFVERSRAYSNCYPYDVPDYVSLRSLVASSGTLEFKNESFNWTGVKQNGTSSACIRRSSSSFFSRLNWLTHLNYTYPALNVTMPNNEQFDKLYIWGIHHPGTDKDQIFLYAQSSGRITVSTKRSQQAVIPNIGSRPRIRDIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKIRSGKSSIMKSDAPIGKCKSECITPNGSISNDKPFQNVNRITYGACPRYVKHSTLKLATGMRNIPEKQTRGIFGAIAGFIENGWEGMKDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDSEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACIGSIRNETYDHNVYRDEALNNRFQIKGIELKSGYKDWILWISFAISCFLLCIALLGFIMWACQKGNIKCNICI",
"MKTIIALSHILCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDRIEVTNATELVQSSSTGEICDSPHHILDGENCTLIDALLGDPQCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLIASSGTLEFNNESFNWAGVTQNGTSSSCIRGSNSSFFSRLNWLTHLNSKYPALNVTMPNKEQFDKLYIWGVHHPGTDKDQISLYAQSSGRITVSTKRSQQAVTPNIGSRPRIRNIPSRISIYWTIVRPGDILLINSTGNLIAPRGYFKIRSGKSSIMRSDAPIGKCNSECITPNGSISNDKPFQNVNKITYGACPRYIKQSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMMDGWYGFRHQNSEGRGQAADLKSTQAAINQITGKLNRVIKKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIDLTDAEMNKLFEKTKKQLRENAEDMGNGCFKIYHKCDNACMGSIRNGTYDHNVYRDEALNNRFQIKGVELKSGYKDWILWISFATSCFLLCVALLGFIMWACQKGNIRCNICI",
"MKTIIALSCILCLVFAQELPGNDDSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGEICNSPHQILDGENCTLMDALLGDPQCDGFQNNKWDLFVERSKAHSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSNNSFFSRLNWLTHLNFKYPALNVTMPNKEQFDKLYIWGVHHPGTDNDQIFLYAQASGRITVSTKRSQQTVIPNIGSRPRVRNIPSRVSIYWTIVKPGDILLINSTGNLIAPRGYFKIRNGKSSIMRSDAPIGKCNSACITPNGSIPNDKPFQNVNRITYGACPRYVKQSTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGLVDGWYGFRHQNSEGTGQAADLKSTQAAIDQINGKLNRLIGKTNEKFHQIEKEFSEVEGRVQDLEKYVEDTKIDLWSYNAELLVALENQHTIHLTDSEMNKLFERTRRQLRENAEDMGNGCFTIYHKCDNACIGSIRNGTYDHDVYRDEALNNRFQIKGVELKSEYKDWILWISFAISCFLLCVALLGFIMWACQKGNIRCNICI"]


def cal_coverage1(to_cal, test_set):
    epitope_nature = [seq[i:i+epitope_length] for seq in test_set for i in range(len(seq)-epitope_length)]
    epitope_mosaic = [seq[i:i+epitope_length] for seq in to_cal for i in range(len(seq)-epitope_length)]
    covers = [epi for epi in epitope_nature if epi in epitope_mosaic]
    return len(covers) / len(epitope_nature)

# averaged for sequence in test set

def cal_coverage2(to_cal, test_set):
    epitope_mosaic = set([seq[i:i+epitope_length] for seq in to_cal
                                                  for i in range(len(seq)-epitope_length)])
    epitope_nature = [[seq[i:i+epitope_length] for i in range(len(seq)-epitope_length)] for seq in test_set]
    covers = []
    for seq in epitope_nature:
        cover = [epi for epi in seq if epi in epitope_mosaic]
        covers.append(len(cover)/len(seq))
    return sum(covers)/len(covers)

# averaged for sequence in mosaic
def cal_coverage3(to_cal, test_set):
    epitope_nature = [seq[i:i+epitope_length] for seq in test_set for i in range(len(seq)-epitope_length)]
    epitope_mosaics = [[seq[i:i+epitope_length] for i in range(len(seq)-epitope_length)] for seq in to_cal]
    covers = []
    for epitope_mosaic in epitope_mosaics:
        single_cover = len([epi for epi in epitope_nature if epi in epitope_mosaic])/len(epitope_nature)
        covers.append(single_cover)
    return mean(covers)

# nr
def cal_coverage4(to_cal, test_set):
    epitope_nature = set([seq[i:i+epitope_length] for seq in test_set for i in range(len(seq)-epitope_length)])
    epitope_mosaic = set([seq[i:i+epitope_length] for seq in to_cal for i in range(len(seq)-epitope_length)])
    covers = [epi for epi in epitope_nature if epi in epitope_mosaic]
    return len(covers) / len(epitope_nature)


def cal_coverage5(to_cal, test_set):
    epitope_mosaic = set([seq[i:i+epitope_length] for seq in to_cal
                                                  for i in range(len(seq)-epitope_length)])
    epitope_nature = [[seq[i:i+epitope_length] for i in range(len(seq)-epitope_length)] for seq in test_set]
    covers = []
    for seq in epitope_nature:
        cover = [epi for epi in seq if epi in epitope_mosaic]
        covers.append(len(cover)/len(seq))
    return sum(covers)/len(covers)


# with open("/home/zeng/Desktop/back.seq", "w") as f:
#     for i, item in enumerate(seq_list):
#         f.write(">back{}".format(i) + "\n" + item + "\n")

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
epitope_nature = [seq[i:i+epitope_length] for seq in seq_list
                                            for i in range(len(seq)-epitope_length)]
test_set = [[seq[i:i+epitope_length] for i in range(len(seq)-epitope_length)]
                                        for seq in seq_list]

epitope_count = Counter(epitope_nature)
no_raw_epitope = {k for k, v in epitope_count.items() if v > 3}   # 超参数 raw number


# print(cal_coverage1(mosaic2, seq_list))
print(cal_coverage2(mosaic2, seq_list))
# # print(cal_coverage3(mosaic2, seq_list))
# # print(cal_coverage4(mosaic2, seq_list))
# print(cal_coverage5(mosaic2, seq_list))


# with open("/home/zeng/Desktop/all_n2.fasta", "w") as f:
#     for i, seq in enumerate(seq_list):
#         f.write(">seq{}\n{}\n".format(i, seq))
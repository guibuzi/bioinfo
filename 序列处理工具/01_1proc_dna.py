import pandas as pd


def import_seq(file_name):
    data = []

    with open(file_name) as f:
        contents = f.read()

    seq_lists = contents.split(sep='>')[1:]
    for seq_list in seq_lists:
        item = {'id': '', 'protein': '', 'name': '', 'epi': '', 'seq': ''}   # 只提取 id 蛋白 及 序列
        a = seq_list.split('\n')
        attr = [x.strip() for x in a[0].split('|')]
        sep = '\n'
        seq = sep.join(a[1:])

        item['id'] = attr[0]
        item['protein'] = attr[1]
        item['name'] = attr[2]
        item['epi'] = attr[3]
        item['seq'] = seq
        data.append(item)

    return data


a = []
for i in range(9):
    tmp = []
    tmp = import_seq(r"暑期实习\\H3N2_to2019_global\\original data\\dna\\dna_file ({0}).fasta".format(i+1))
    a.append(tmp)

# 合并
dup_total = []

for i in a:
    dup_total.extend(i)


# 转化为pandas结构
frame = pd.DataFrame(dup_total)

# 去重
frame2 = frame.drop_duplicates()

# 统计蛋白种类
unique_protein = frame2['protein'].unique()


def export_seq(pro):
    frame3 = frame2[frame2['protein'] == pro]
    frame4 = frame3.drop_duplicates('epi')
    with open(r"D:\\python_work\\dna\\dna_{0}.fasta".format(pro), 'w') as f2:
        for i in frame4.index:
            name = frame4.at[i, 'id']
            seq = frame4.at[i, 'seq']
            f2.write(">EPI" + name + "\n" + seq)


for i in unique_protein:
    export_seq('{}'.format(i))

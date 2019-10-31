# 读入fasta文件
data = []

with open('protein2004.fasta') as f:
    contents = f.read()

seq_lists = contents.split(sep='>')[1:]
for seq_list in seq_lists:
    item = {'id': '', 'protein': '', 'name': '', 'epi': '', 'loci': '', 'type': '', 'seq': ''}
    a = seq_list.split('\n')
    attr = [x.strip() for x in a[0].split('|')]
    sep = ''
    seq = sep.join(a[1:])

    item['id'] = attr[0]
    item['protein'] = attr[1]
    item['name'] = attr[2]
    item['epi'] = attr[3]
    item['loci'] = attr[4]
    item['type'] = attr[5]
    item['seq'] = seq
    data.append(item)


# 计算不同蛋白类型
protein = []
for seq_list in data:
    a = seq_list['protein']
    protein.append(a)

unique_protein = set(protein)


# 按蛋白类型写出文件

for protein in unique_protein:
    with open("protein_{0}_2004.fasta".format(protein), 'w') as f2:
        for seq_list in data:
            name = seq_list['id']
            seq = seq_list['seq']
            f2.write(">" + name + "\n" + seq + '\n')

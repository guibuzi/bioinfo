def import_seq(file_name):
    data = []

    with open(file_name) as f:
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

    return data


# 读入文件
a = import_seq(file_name='protein2004.fasta')
a1 = import_seq(file_name='protein.fasta')


# 合并去重
combination_data = a[0:10] + a1
unique_combination = []

for seq_list in combination_data:
    if seq_list not in unique_combination:
        unique_combination.append(seq_list)

print(unique_combination)

from collections import Counter, defaultdict

# 读入fasta文件
def read_fasta(file_path):
    """
    describtion: 读取 fasta 文件

    parameters:
    file_path: fasta 文件所在路径
    
    return: dict of seq
    """
    data = []
    with open(file_path) as f:
        contents = f.read()

    seq_lists = contents.split(sep='>')[1:]
    for seq_list in seq_lists:
        item = defaultdict(str)

        a = seq_list.split('\n')
        attr = [x.strip() for x in a[0].split('|')]
        seq = "".join(a[1:]) 
        item['id'] = attr[0]
        item['protein'] = attr[1]
        item['name'] = attr[2]
        item['epi'] = attr[3]
        item['seq'] = seq
        data.append(item)
    return data


def split_by_segment(data):
    # 计算不同蛋白类型
    protein = []
    uni_pro = []
    for seq_list in data:
        a = seq_list['protein']
        protein.append(a)

    counter_by_segment = Counter(protein)
    print("\n============================== 本轮更新的序列有 ==============================")
    for key, value in counter_by_segment.items():
        uni_pro.append(key)
        print("The number of {0} is: {1}".format(key, str(value)))
    print("================================ 本轮更新完毕 ================================\n")
    
    # 按蛋白类型写出文件
    for protein in uni_pro:
        with open("/home/zeng/Desktop/H3N2/data/Protein/protein_{}.fasta".format(protein), 'a') as fw:
            for seq_list in data:
                if protein == seq_list['protein']:
                    name = seq_list['id']
                    epi = seq_list['epi']
                    seq = seq_list['seq']
                    fw.write(">" + name + "\n" + seq + '\n')


def update_fasta(data_path, old_entry):
    data = read_fasta(data_path)
    new_data = [x for x in data if x['epi'] not in old_entry]
    new_entry = set(x['epi'] for x in data if x['epi'] not in old_entry)
    old_entry.update(new_entry)
    split_by_segment(new_data)
    return old_entry


if __name__ == "__main__":
    old_entry = set()
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (1).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (2).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (3).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (4).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (5).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (6).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (7).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (8).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/Desktop/H3N2/data/original_proteion/protein_file (9).fasta", old_entry)

    current_entry = old_entry

    current_entry = update_fasta("/home/zeng/Desktop/US_oct/gisaid_epiflu_sequence.fasta", current_entry)

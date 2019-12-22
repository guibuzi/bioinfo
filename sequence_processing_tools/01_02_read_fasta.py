from collections import Counter, defaultdict
import json


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
        item['seq'] = seq
        data.append(item)
    return data


def split_by_segment(data):
    # 计算不同蛋白类型
    uni_pro = []
    protein =[seq_list['protein'] for seq_list in data]

    counter_by_segment = Counter(protein)
    print("============================== 本轮更新的序列有 ==============================")
    for key, value in counter_by_segment.items():
        uni_pro.append(key)
        print("The number of {0} is: {1}".format(key, str(value)))
    print("================================ 本轮更新完毕 ================================\n")

    # 按蛋白类型写出文件
    for protein in uni_pro:
        with open("/home/zeng/INFLUENZA_DATABASE/H3N2/data/Protein/protein_{}.fasta".format(protein), 'a') as fw:
            for seq_list in data:
                if protein == seq_list['protein']:
                    name = seq_list['id']
                    seq = seq_list['seq']
                    fw.write(">" + name + "\n" + seq + '\n')


def update_fasta(data_path, old_entry):
    print("There are already have %s segments sequeces." % len(old_entry))
    data = read_fasta(data_path)
    new_data = [x for x in data if x['id'] not in old_entry]
    new_entry = set(x['id'] for x in data if x['id'] not in old_entry)
    print("%s new segments sequences will be added." % len(new_entry))
    old_entry.update(new_entry)       # renew the set to restore the segment epi
    split_by_segment(new_data)
    return old_entry


def read_old(isl_path):
    with open(isl_path) as f:
        contents = f.readlines()
    if contents:
        return set(contents)
    else:
        return set()


if __name__ == "__main__":
    # read the old json file
    old_entry = read_old("/home/zeng/INFLUENZA_DATABASE/H3N2/data/Protein/segment_id.txt")

    # change this line when add new isolation
    # you can add new segment sequences from multiple fasta file
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (1).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (2).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (3).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (4).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (5).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (6).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (7).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (8).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/protein_file (9).fasta", old_entry)
    old_entry = update_fasta("/home/zeng/INFLUENZA_DATABASE/H3N2/data/original/Protein/gisaid_epiflu_sequence.fasta", old_entry)

    # when finish read all the new isolations, re-write the json file
    with open("/home/zeng/INFLUENZA_DATABASE/H3N2/data/Protein/segment_id.txt", "w") as f:
        f.writelines(list(old_entry))

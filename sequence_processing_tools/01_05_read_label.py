import os, sys
import json


def read_label(path):
    with open(path) as f:
        lines = f.readlines()
    return {line.split()[0].strip(): line.split()[1].strip() for line in lines}

def read_info(path):
    with open(path) as f:
        data = json.loads(f.read())
    return data

def segment_to_epi(data):
    return {value['HA']:key for key, value in data.items()}

def read_segment(path):
    data = {}
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                label = line[1:-1]
                data[label] = ''
            else:
                data[label] += line.replace("\n", "")
    return data

def read_index(path):
    with open(path) as f:
        data = f.readlines()
    return [i.replace("\n", "") for i in data]


if __name__ == "__main__":
    label_path = "/home/zeng/Desktop/Global_traced_label/"
    data_path = "/home/zeng/INFLUENZA_DATABASE/H3N2/data/"
    sequence_path = data_path + "Protein/"

    # all
    all_genomic = read_index("/home/zeng/Desktop/index")

    # read isolation information
    info_data  =read_info(data_path + "isolation_information.json")
    segment_2_epi = segment_to_epi(info_data)

    # read sequences of all segments
    all_gene = {}
    protein_file =[item for item in os.listdir(sequence_path) if item.endswith(".fasta")]
    segments = [item.split("_")[1].replace(".fasta", "") for item in protein_file]
    
    for _seqment, _file in zip(segments, protein_file):
        current_path = os.path.join(sequence_path, _file)
        data = read_segment(current_path)
        all_gene[_seqment] = data

    # read ha epi with cluster information from all regions
    global_label = {}
    files = os.listdir(label_path)
    regions = [name.split("_")[0] for name in files]

    for _region, _file in zip(regions, files):
        current_path = os.path.join(label_path, _file)
        data = read_label(current_path)
        global_label[_region] = data
    
    # search the sequence for each segment each region
    out_path = "/home/zeng/Desktop/reassortment_analysis/"

    not_found = []
    not_found2 = []
    not_found3 = []
    not_found4 = []

    for _region in regions:
        region_path = out_path + _region
        try:
            os.mkdir(region_path)
        except:
            pass

        # read isl id from ha id
        labels = global_label[_region]

        isl_labels = []
        for _label in labels.keys():
            try:
                if segment_2_epi[_label] in all_genomic:
                    isl_labels.append(segment_2_epi[_label])
                else:
                    print("Drop %s in %s due to lack segments." % (_label, _region))
                    not_found4.append(" ".join([_region, _label, "\n"]))
            except:
                # 本地数据库不全 输出缺少的ha ===》 去gisaid找
                print("DATABASE Error: Not found HA %s in %s !" % (_label, _region))
                not_found.append(" ".join([_region, _label, '\n']))
        

        segment_map = {"PB2": "PB2", "PB1": "PB1", "PA": "PA", "HA": "HA", "NP":"NP", "NA": "NA", "M1": "MP", "M2": "MP", "NS1":"NS", "NS2": "NS", "PB1-F2": "PB1", "NEP": "NS"}
        
        for _isl_label in isl_labels:
            date = info_data[_isl_label]['date']
            name = info_data[_isl_label]['name']
            cluster = labels[info_data[_isl_label]['HA']]
            for _segment, _segment_key in segment_map.items():
                _segment_epi = info_data[_isl_label][_segment_key]
                if _segment_epi:
                    try:
                        sequence = all_gene[_segment][_segment_epi]
                    except:
                        # 本地库不全，有片段id 但是本地库缺少对应序列
                        if _segment != "NS2":
                            print("DATABASE Error: Not fount %s 's gene %s (%s) in %s" % (_isl_label, _segment, _segment_epi, _region))
                            not_found3.append(" ".join([_isl_label, _segment, _segment_epi, _region, "\n"]))
                else:
                    # 网络库不全，gisaid 未包括的片段
                    print("SEGMENTS Error: Not fount %s 's gene %s in %s" % (_isl_label, _segment, _region))
                    not_found2.append(" ".join([_isl_label, _segment_key, _region, "\n"]))

                with open(os.path.join(region_path, "%s.fasta" % _segment), "a") as g:
                    g.write(">%s|%s|%s|%s|%s\n%s\n" % (_isl_label, _segment_epi, cluster, name, date, sequence))


    with open("/home/zeng/Desktop/log1", "w") as l:
        l.writelines(not_found)
    
    with open("/home/zeng/Desktop/log2", "w") as l:
        l.writelines(not_found2)

    with open("/home/zeng/Desktop/log3", "w") as l:
        l.writelines(not_found3)

    with open("/home/zeng/Desktop/log4", "w") as l:
        l.writelines(not_found4)

import pandas as pd
import re


def proc_tree(infile, outfile, segment):
    df2 = df.loc[(df['{} Segment_Id'.format(segment)].notna()), ['Isolate_Id',
                                                                 '{} Segment_Id'.format(segment), 'Isolate_Name', 'year']]
    df2['{} Segment_Id'.format(segment)] = df2['{} Segment_Id'.format(segment)].str.strip()
    df2 = df2.set_index(['{} Segment_Id'.format(segment)])

    with open(infile, 'r') as fr:
        _contents = fr.read()
        _regex = re.compile(r'EPI[0-9]*')
        _segment_epi_list = re.findall(_regex, _contents)
        _dict_of_epi = df2[df2.index.isin(_segment_epi_list)].to_dict()

        for epi in _segment_epi_list:
            _regex2 = re.compile(epi + ":")
            _repl = '"' + _dict_of_epi['Isolate_Id'][epi] + '_' + str(_dict_of_epi['year'][epi]) + \
                    '"[&year=' + str(_dict_of_epi['year'][epi]) + ', name="' + _dict_of_epi['Isolate_Name'][epi] + '"]:'
            _contents = re.sub(_regex2, _repl, _contents)

    with open(outfile, 'w', encoding='utf-8') as fw:
        fw.write("#NEXUS\n")
        fw.write("Begin trees;\n")
        fw.write("\ttree Tree1= [&R] " + _contents + "\n")
        fw.write("end;")


if __name__ == "__main__":
    df = pd.read_excel(r"D:\ioslate_information_处理.xlsx")
    df['region'] = df['Location'].str.split(' / ').str.get(0)
    df['country'] = df['Location'].str.split(' / ').str.get(1)
    df['year'] = df['Collection_Date'].str.split('-').str.get(0)
    df = df[(df['region'] == 'Europe')]

    # proc_tree(r"D:\usa\tree\usa_HA_tree", r"D:\usa\tree\usa_HA_tree_re", "HA")
    # proc_tree(r"D:\usa\tree\usa_MP_tree", r"D:\usa\tree\usa_M1_tree_re", "MP")
    # proc_tree(r"D:\usa\tree\usa_M2_tree", r"D:\usa\tree\usa_M2_tree_re", "MP")
    # proc_tree(r"D:\usa\tree\usa_NA_tree", r"D:\usa\tree\usa_NA_tree_re", "NA")
    # proc_tree(r"D:\usa\tree\usa_NP_tree", r"D:\usa\tree\usa_NP_tree_re", "NP")
    # proc_tree(r"D:\usa\tree\usa_NS_tree", r"D:\usa\tree\usa_NS_tree_re", "NS")
    # proc_tree(r"D:\usa\tree\usa_PA_tree", r"D:\usa\tree\usa_PA_tree_re", "PA")
    # proc_tree(r"D:\usa\tree\usa_PB1_tree", r"D:\usa\tree\usa_PB1_tree_re", "PB1")
    # proc_tree(r"D:\usa\tree\usa_PB2_tree", r"D:\usa\tree\usa_PB2_tree_re", "PB2")

    # proc_tree(r"D:\australia\tree\Australia_HA_tree", r"D:\australia\tree\Australia_HA_tree_re", "HA")
    # proc_tree(r"D:\australia\tree\Australia_M1_tree", r"D:\australia\tree\Australia_M1_tree_re", "MP")
    # proc_tree(r"D:\australia\tree\Australia_M2_tree", r"D:\australia\tree\Australia_M2_tree_re", "MP")
    # proc_tree(r"D:\australia\tree\Australia_NA_tree", r"D:\australia\tree\Australia_NA_tree_re", "NA")
    # proc_tree(r"D:\australia\tree\Australia_NP_tree", r"D:\australia\tree\Australia_NP_tree_re", "NP")
    # proc_tree(r"D:\australia\tree\Australia_NS_tree", r"D:\australia\tree\Australia_NS_tree_re", "NS")
    # proc_tree(r"D:\australia\tree\Australia_PA_tree", r"D:\australia\tree\Australia_PA_tree_re", "PA")
    # proc_tree(r"D:\australia\tree\Australia_PB1_tree", r"D:\australia\tree\Australia_PB1_tree_re", "PB1")
    # proc_tree(r"D:\australia\tree\Australia_PB2_tree", r"D:\australia\tree\Australia_PB2_tree_re", "PB2")

    # proc_tree(r"D:\china\tree\China_HA_tree", r"D:\china\tree\China_HA_tree_re", "HA")
    # proc_tree(r"D:\china\tree\China_M1_tree", r"D:\china\tree\China_M1_tree_re", "MP")
    # proc_tree(r"D:\china\tree\China_M2_tree", r"D:\china\tree\China_M2_tree_re", "MP")
    # proc_tree(r"D:\china\tree\China_NA_tree", r"D:\china\tree\China_NA_tree_re", "NA")
    # proc_tree(r"D:\china\tree\China_NP_tree", r"D:\china\tree\China_NP_tree_re", "NP")
    # proc_tree(r"D:\china\tree\China_NS_tree", r"D:\china\tree\China_NS_tree_re", "NS")
    # proc_tree(r"D:\china\tree\China_PA_tree", r"D:\china\tree\China_PA_tree_re", "PA")
    # proc_tree(r"D:\china\tree\China_PB1_tree", r"D:\china\tree\China_PB1_tree_re", "PB1")
    # proc_tree(r"D:\china\tree\China_PB2_tree", r"D:\china\tree\China_PB2_tree_re", "PB2")

    # proc_tree(r"D:\europe\tree\Europe_HA_tree", r"D:\europe\tree\Europe_HA_tree_re", "HA")
    # proc_tree(r"D:\europe\tree\Europe_M1_tree", r"D:\europe\tree\Europe_M1_tree_re", "MP")
    # proc_tree(r"D:\europe\tree\Europe_M2_tree", r"D:\europe\tree\Europe_M2_tree_re", "MP")
    # proc_tree(r"D:\europe\tree\Europe_NA_tree", r"D:\europe\tree\Europe_NA_tree_re", "NA")
    # proc_tree(r"D:\europe\tree\Europe_NP_tree", r"D:\europe\tree\Europe_NP_tree_re", "NP")
    # proc_tree(r"D:\europe\tree\Europe_NS_tree", r"D:\europe\tree\Europe_NS_tree_re", "NS")
    # proc_tree(r"D:\europe\tree\Europe_PA_tree", r"D:\europe\tree\Europe_PA_tree_re", "PA")
    # proc_tree(r"D:\europe\tree\Europe_PB1_tree", r"D:\europe\tree\Europe_PB1_tree_re", "PB1")
    # proc_tree(r"D:\europe\tree\Europe_PB2_tree", r"D:\europe\tree\Europe_PB2_tree_re", "PB2")

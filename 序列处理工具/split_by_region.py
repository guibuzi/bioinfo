import pandas as pd


def substrate_region(region, df):
    _df2 = df[df['region'].isin([region])]

    _protein_list = {'HA': 'HA', 'NA': 'NA', 'PA': 'PA', 'PB1': 'PB1', 'PB2': 'PB2',
                     'M1': 'MP', 'M2': 'MP', 'NP': 'NP', 'NS': 'NS'}
    for _key, _value in _protein_list.items():
        _df3 = _df2[_df2['{} Segment_Id'.format(_value)].notna()][['Isolate_Id',
                                                                   '{} Segment_Id'.format(_value), 'Isolate_Name', 'year']]
        _df3['{} Segment_Id'.format(_value)] = _df3['{} Segment_Id'.format(_value)].str.strip()
        segment_list = _df3['{} Segment_Id'.format(_value)].to_list()

        _protein = r'D:\protein_{}.fasta'.format(_key)
        _sub_protein = r'D:\{0}_{1}.fasta'.format(region, _key)

        with open(_protein, 'r') as f:
            _segment_dict = {}
            _lines = f.readlines()
            _ha_epi = ''
            for _line in _lines:
                if _line.startswith('>'):
                    _ha_epi = _line[1:-1]
                else:
                    _segment_dict[_ha_epi] = _line[:-1]

        _want = {ha: _segment_dict.get(ha) for ha in segment_list if _segment_dict.get(ha)}

        with open(_sub_protein, 'w') as fw:
            for key, value in _want.items():
                fw.write(">" + key + "\n")
                fw.write(value + "\n")


if __name__ == '__main__':
    df = pd.read_excel(r"D:\python_work\gisadi_data\ioslate_information_处理.xlsx")
    df['region'] = df['Location'].str.split(' / ').str.get(0)
    df['country'] = df['Location'].str.split(' / ').str.get(1)
    df['year'] = df['Collection_Date'].str.split('-').str.get(0)
    substrate_region('Europe', df)

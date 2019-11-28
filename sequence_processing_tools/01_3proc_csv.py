import pandas as pd


def import_xls(i):
    df = pd.read_excel(r"暑期实习\\csv\\ioslate_information ({0}).xls".format(i))

    df['PB2 Segment_Id'] = df['PB2 Segment_Id'].str.split('|').str.get(0)
    df['PB1 Segment_Id'] = df['PB1 Segment_Id'].str.split('|').str.get(0)
    df['PA Segment_Id'] = df['PA Segment_Id'].str.split('|').str.get(0)
    df['HA Segment_Id'] = df['HA Segment_Id'].str.split('|').str.get(0)
    df['NP Segment_Id'] = df['NP Segment_Id'].str.split('|').str.get(0)
    df['NA Segment_Id'] = df['NA Segment_Id'].str.split('|').str.get(0)
    df['MP Segment_Id'] = df['MP Segment_Id'].str.split('|').str.get(0)
    df['NS Segment_Id'] = df['NS Segment_Id'].str.split('|').str.get(0)
    # df['HE Segment_Id'] = df['HE Segment_Id'].str.split('|').str.get(0)
    # df['P3 Segment_Id'] = df['P3 Segment_Id'].str.split('|').str.get(0)

    return df


df1 = import_xls(1)
df2 = import_xls(2)
df3 = import_xls(3)
df4 = import_xls(4)
df5 = import_xls(5)
df6 = import_xls(6)
df7 = import_xls(7)
df8 = import_xls(8)
df9 = import_xls(9)

frame = [df1, df2, df3, df4, df5, df6, df7, df8, df9]
result = pd.concat(frame)
unique_result = result.drop_duplicates()

unique_result.to_excel(r'暑期实习\\csv\\ioslate_information_处理.xlsx', sheet_name='Sheet1')

def to_one_row(infile, outfile):
    fr = open(infile, 'r')   # 输入文件
    fw = open(outfile, 'w')  # 输出文件
    for line in fr:
        if line.startswith('>'):    # 判断字符串是否以‘>开始’
            fw.write("\n" + line)
        else:
            line = line.replace('\n', '')
            fw.write(line)
    fr.close()
    fw.close()


if __name__ == '__main__':
    to_one_row(r"D:\europe\alignment\Europe_HA_re.fasta",
               r"D:\europe\alignment\Europe_HA_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_M1_re.fasta",
               r"D:\europe\alignment\Europe_M1_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_M2_re.fasta",
               r"D:\europe\alignment\Europe_M2_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_NA_re.fasta",
               r"D:\europe\alignment\Europe_NA_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_NP_re.fasta",
               r"D:\europe\alignment\Europe_NP_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_NS_re.fasta",
               r"D:\europe\alignment\Europe_NS_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_PA_re.fasta",
               r"D:\europe\alignment\Europe_PA_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_PB1_re.fasta",
               r"D:\europe\alignment\Europe_PB1_re2.fasta")
    to_one_row(r"D:\europe\alignment\Europe_PB2_re.fasta",
               r"D:\europe\alignment\Europe_PB2_re2.fasta")

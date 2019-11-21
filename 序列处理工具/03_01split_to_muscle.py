def split_file(to_split, out_path, limit=2000):
    file_count = 0
    url_list = []

    with open(to_split) as f:
        for line in f:
            url_list.append(line)
            if len(url_list) < limit:
                continue
            # 数据达到LIMIT
            file_name = str(file_count)+ ".fasta"
            with open(out_path + file_name, 'w') as file:
                for url in url_list[:-1]:
                    #print(url)
                    file.write(url)
                file.write(url_list[-1].strip())
            url_list = []
            file_count += 1
        if url_list:
            file_name=str(file_count) + ".fasta"
            with open(out_path + file_name, 'w') as file:
                for url in url_list:
                    file.write(url)


split_file("/home/zeng/Desktop/US_oct/data/usa_HA.fasta", "/home/zeng/Desktop/US_oct/data/ha/")
split_file("/home/zeng/Desktop/US_oct/data/usa_M1.fasta", "/home/zeng/Desktop/US_oct/data/m1/")
split_file("/home/zeng/Desktop/US_oct/data/usa_NA.fasta", "/home/zeng/Desktop/US_oct/data/na/")
split_file("/home/zeng/Desktop/US_oct/data/usa_NP.fasta", "/home/zeng/Desktop/US_oct/data/np/")
split_file("/home/zeng/Desktop/US_oct/data/usa_NS1.fasta", "/home/zeng/Desktop/US_oct/data/ns1/")
split_file("/home/zeng/Desktop/US_oct/data/usa_PA.fasta", "/home/zeng/Desktop/US_oct/data/pa/")
split_file("/home/zeng/Desktop/US_oct/data/usa_PB1.fasta", "/home/zeng/Desktop/US_oct/data/pb1/")
split_file("/home/zeng/Desktop/US_oct/data/usa_PB2.fasta", "/home/zeng/Desktop/US_oct/data/pb2/")

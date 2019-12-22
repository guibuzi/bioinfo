def split_file(to_split, out_path, limit=3000):
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


split_file("/home/zeng/Desktop/no_redundant_h3.fasta", "/home/zeng/Desktop/no_redundant_h3/")

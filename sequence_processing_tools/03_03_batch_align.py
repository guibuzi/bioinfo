import sys, os
import multiprocessing


def split_file(to_split, out_split, limit=4000):
    file_count = 0
    url_list = []
    
    with open(to_split) as f:
        for line in f:
            url_list.append(line)
            if len(url_list) < limit:
                continue
            # æ•°æ�®è¾¾åˆ°LIMIT
            file_name = str(file_count)+ ".fasta"
            with open(out_split + "/" + file_name, 'w') as file:
                for url in url_list[:-1]:
                    #print(url)
                    file.write(url)
                file.write(url_list[-1].strip())
            url_list = []
            file_count += 1
        if url_list:
            file_name=str(file_count) + ".fasta"
            with open(out_split + "/" + file_name, 'w') as file:
                for url in url_list:
                    file.write(url)


if __name__ == "__main__":
    # change this part to run
    tasks = ['Aus', 'Europe', 'SEAsia', 'USA']
    path = "/home/zeng/Desktop/reassortment_analysis/"

    for task in tasks:
        try:
            os.makedirs(path + "out_split/" + task)
            os.makedirs(path + "out_align/" + task)
        except:
            os.system("rm -rf %s" % (path + "out_split/"))
            os.system("rm -rf %s" % (path + "out_align/"))
            os.makedirs(path + "out_split/" + task)
            os.makedirs(path + "out_align/" + task)

    for task in tasks:
        to_align_dir = path + task
        to_align = os.listdir(to_align_dir)
        for i in to_align:
            file_name = i.split(".")[0]
            os.mkdir(path + "out_align/" + task + "/" + file_name)
            os.mkdir(path + "out_split/" + task + "/" + file_name)
            split_file(path + task + "/" + i, path + "out_split/" + task + "/" + file_name)


    in_seqs = []
    out_seqs = []

    for task in tasks:
        to_align_dir = path + "out_split/" + task
        to_path = path + "out_align/" + task
        to_align_segment = os.listdir(to_align_dir)
        for segment in to_align_segment:
            to_align = os.listdir(to_align_dir + "/" + segment)
            for j in to_align:
                in_seqs.append(os.path.join(to_align_dir + "/" + segment, j))
                out_seqs.append(os.path.join(to_path + "/" + segment, j))

    for i, o in zip(in_seqs, out_seqs):
        print(i, o)


    def muscle(_in, _out):
        print("Begin to align %s" % os.path.split(_in)[1])
        os.system("muscle -in %s -out %s -quiet" % (_in, _out))
        print("Finsh align %s" % os.path.split(_in)[1])


    pool = multiprocessing.Pool(processes = 6)
    print("Begin tasks ...")
    for i, o in zip(in_seqs, out_seqs):
        pool.apply_async(muscle, (i, o))

    pool.close()
    pool.join()
    print("Finsh muscle.")

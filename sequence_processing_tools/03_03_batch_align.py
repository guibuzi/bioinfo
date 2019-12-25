import sys, os
import multiprocessing


def split_file(to_split, out_path, limit=4000):
    file_count = 0
    url_list = []

    with open(to_split) as f:
        for line in f:
            url_list.append(line)
            if len(url_list) < limit:
                continue
            # æ•°æ�®è¾¾åˆ°LIMIT
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




# change this part to run
tasks = ['Aus', 'Europe', 'SEAsia']
current_path = "/home/zeng/Desktop/reassortment_analysis/"

for task in tasks:
    try:
        os.makedirs(current_path + "align/" + task)
    except:
        os.system("rm -rf %s" % (current_path + "align/"))
        os.makedirs(current_path + "align/" + task)

in_seqs = []
out_seqs = []

for task in tasks:
    path = current_path + task
    to_path = current_path + "align/" + task
    to_align = os.listdir(path)
    for f in to_align:
        in_seqs.append(os.path.join(path, f))
        out_seqs.append(os.path.join(to_path, f))

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

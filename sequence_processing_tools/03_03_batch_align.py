import sys, os
import multiprocessing


tasks = ['BY', 'BV', 'H1-A', 'H3']
current_path = "/home/zeng/Desktop/cao/"

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


pool = multiprocessing.Pool(processes = 4)
print("Begin tasks ...")
for i, o in zip(in_seqs, out_seqs):
    pool.apply_async(muscle, (i, o))

pool.close()
pool.join()
print("Finsh muscle.")

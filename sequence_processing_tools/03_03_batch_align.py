import sys, os
from multiprocessing import Pool


def muscle(_in, _out):
    print("Begin to align %s" % os.path.split(_in)[1])
    os.system("muscle -in %s -out %s -quiet" % (_in, _out))
    print("Finsh align %s" % os.path.split(_in)[1])


path = "/home/zeng/Desktop/sample2"
out_path = "/home/zeng/Desktop/sample2/align"

os.system("rm -rf /home/zeng/Desktop/sample2/align")

tasks = os.listdir(path)
in_seqs = []
out_seqs = []

for task in tasks:
    _in = os.path.join(path,  task)
    _out = os.path.join(out_path, task)
    in_seqs.append(_in)
    out_seqs.append(_out)

os.makedirs("/home/zeng/Desktop/sample2/align", exist_ok=True)

pool = Pool(processes = 6)
print("Begin tasks ...")
for i, o in zip(in_seqs, out_seqs):
    pool.apply_async(muscle, (i, o))

pool.close()
pool.join()
print("Finsh muscle.")

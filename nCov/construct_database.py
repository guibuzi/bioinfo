import sys, os

def blastdb_aliastool(_gi, _out):
    command = "blastdb_aliastool -db nt/nt -gilist %s -dbtype nucl -out %s" % (_gi, _out)
    os.system(command)

def blastdbcmd(_in, _out):
    command = "blastdbcmd -db %s -entry all -out %s" % (_in, _out)
    os.system(command)


os.chdir("/home/zeng/blastdb/coronaviridae")

path = "/home/zeng/Desktop/gi_dir"
tasks = os.listdir(path)
in_list = [os.path.join(path, task) for task in tasks]
db_list = [task.split('.')[0] for task in tasks]
fasta_list = [db + ".fasta" for db in db_list]

for _in, _out in zip(in_list, db_list):
    print(_in, _out)
    # blastdb_aliastool(_in, _out)

for _in, _out in zip(db_list, fasta_list):
    print(_in, _out)
    blastdbcmd(_in, _out)
#import modules
import os
import argparse
import pandas as pd

#2022/07/06 created by yuanw

#print help information
parser = argparse.ArgumentParser(description='Divide subgenomes.')
parser.add_argument("-sp4gff",help="input allotetraploid species gff file")
parser.add_argument("-sp4pep",help="input allotetraploid species pep file")
parser.add_argument("-sp4cds",help="input allotetraploid species cds file")
parser.add_argument("-sp4genome",help="input allotetraploid genome fasta file")
args = parser.parse_args()

#dir
path = os.getcwd()
path2 = path + "/2.Cluster"
path3 = path + "/3.Tree"
os.mkdir(path + '/4.groupAB')
path4 = path + "/4.groupAB"

#identity subA or subB
os.system("grep '(sp2,sp41)' %s|cut -f 1 > %s" % (path3 + "/all.treefile.bootstrap",path4 + "/group_A_B.txt"))
os.system("grep '(sp2,sp42)' %s|cut -f 1 > %s" % (path3 + "/all.treefile.bootstrap",path4 + "/group_B_A.txt"))
os.system("grep '(sp41,sp42)' %s|cut -f 1 > %s" % (path3 + "/all.treefile.bootstrap",path4 + "/group_NS.txt"))

f1 = open(path4 + "/group_A_B.txt","r")
f2 = open(path4 + "/group_B_A.txt","r")
f3 = open(args.sp4gff,"r")
f4 = open(path2 + "/sp4_sp2.1v2.cluster.number","r")
f5 = open(path4 + "/groupAB.result","w")

lst1 = []
for line1 in f1.readlines():
    lst1.append(line1.strip())

lst2 = []
for line2 in f2.readlines():
    lst2.append(line2.strip())

name_dict = {}
for line in f3.readlines():
    k = line.strip().split('\t')[1]
    lst = []
    lst.append(str(line.strip().split('\t')[0]))
    lst.append(str(line.strip().split('\t')[2]))
    lst.append(str(line.strip().split('\t')[3]))
    name_dict[k] = lst

for line3 in f4.readlines():
    if line3.strip().split("\t")[0] in lst1:
        f5.write("A" + "\t" + line3.strip().split("\t")[2] + "\t" + '\t'.join(name_dict[line3.strip().split("\t")[2]]) + "\n" +
                 "B" + "\t" + line3.strip().split("\t")[3] + "\t" + '\t'.join(name_dict[line3.strip().split("\t")[3]])  + "\n")
    elif line3.strip().split("\t")[0] in lst2:
        f5.write("B" + "\t" + line3.strip().split("\t")[2] + "\t" + '\t'.join(name_dict[line3.strip().split("\t")[2]]) + "\n" +
                 "A" + "\t" + line3.strip().split("\t")[3] + "\t" + '\t'.join(name_dict[line3.strip().split("\t")[3]])  + "\n")
    else:
        continue

f1.close()
f2.close()
f3.close()
f4.close()
f5.close()

#countAB in scaffold
f1 = open(path4 + "/groupAB.result","r")
f2 = open(path4 + "/groupAB.gene.name","w")
for line in f1.readlines():
	line = line.strip().split('\t')
	f2.write(str(line[0]) + "_" + str(line[2]) + "\n")
f1.close()
f2.close()

df = pd.read_csv(path4 + "/groupAB.gene.name",sep = '\t',header = None, names = ['id'])
df0 = df.groupby(['id'])['id'].count()
df0.to_csv(path4 + "/groupAB.gene.name.count",header=False,sep='\t')

os.system("cut -f 3 %s|sort|uniq> %s" % (path4 + "/groupAB.result",path4 + "/scaffold.name"))

f1 = open(path4 + "/groupAB.gene.name.count","r")
f2 = open(path4 + "/scaffold.name","r")
f3 = open(path4 + "/groupAB.summary","w")

count_dict = {}
for line1 in f1.readlines():
    k = line1.strip().split('\t')[0]
    v = line1.strip().split('\t')[1]
    count_dict[k] = v

for line2 in f2.readlines():
    line2=line2.strip()
    nameA = "A_" + str(line2)
    nameB = "B_" + str(line2)

    if nameA in count_dict:
        countA = count_dict[nameA]
    else:
        countA = 0

    if nameB in count_dict:
        countB = count_dict[nameB]
    else:
        countB = 0

    f3.write(str(line2) + "\t" + "A" + "\t" + str(countA)  + "\t" + "B" + "\t" +  str(countB) + "\n")

f1.close()
f2.close()
f3.close()

#groupAB(A/B or B/A >90%)
os.system("seqkit fx2tab -n -i -l %s > %s" % (args.sp4genome,path4 + "/sp4.genome.fa.len"))

f1 = open(path4 + "/sp4.genome.fa.len","r")
f2 = open(path4 + "/groupAB.summary","r")
f3 = open(path4 + "/group.txt","w")

len_dict = {}
for line1 in f1.readlines():
    k = line1.strip().split('\t')[0]
    v = line1.strip().split('\t')[1]
    len_dict[k] = v

for line2 in f2.readlines():
    countA = line2.strip().split('\t')[2]
    countB = line2.strip().split('\t')[4]
    total = int(countA) + int(countB)
    rate1 = int(countA)/int(total)
    rate2 = int(countB)/int(total)
    len = len_dict[line2.strip().split('\t')[0]]
    if int(countA) > int(countB):
        if rate1 > 0.9:
            f3.write(str(line2.strip().split('\t')[0]) + "\t" + "A" + "\t" + str(len) + "\n")
    elif int(countB) > int(countA):
        if rate2 > 0.9:
            f3.write(str(line2.strip().split('\t')[0]) + "\t" + "B" + "\t" + str(len) + "\n")

f1.close()
f2.close()
f3.close()

#divide
f1 = open(path4 + "/group.txt","r")
f2 = open(args.sp4gff,"r")
f3 = open(path4 + "/subA.scaffold.name","w")
f4 = open(path4 + "/subB.scaffold.name","w")
f5 = open(path4 + "/subA.gene.name","w")
f6 = open(path4 + "/subB.gene.name","w")

lst_A = []
lst_B = []
for line1 in f1.readlines():
    line1 = line1.strip()
    if str(line1.split('\t')[1]) == "A":
        lst_A.append(line1.split('\t')[0])
        f3.write(str(line1.split('\t')[0]) + "\n")
    else:
        lst_B.append(line1.split('\t')[0])
        f4.write(str(line1.split('\t')[0]) + "\n")

for line2 in f2.readlines():
    line2.strip()
    if line2.split('\t')[0] in lst_A:
        f5.write(line2.split('\t')[1] + "\n")
    elif line2.split('\t')[0] in lst_B:
        f6.write(line2.split('\t')[1] + "\n")
    else:
        continue

f1.close()
f2.close()
f3.close()
f4.close()
f5.close()
f6.close()

#result
os.mkdir(path + '/5.result')
path5 = path + "/5.result"

#gff
f1 = open(path4 + "/subA.scaffold.name","r")
f2 = open(path4 + "/subB.scaffold.name","r")
f3 = open(args.sp4gff,"r")
f4 = open(path5 + "/subA.gff","w")
f5 = open(path5 + "/subB.gff","w")

subA_list = []
for line1 in f1.readlines():
	subA_list.append(line1.strip())

subB_list = []
for line2 in f2.readlines():
	subB_list.append(line2.strip())

for line3 in f3.readlines():
	if line3.strip().split('\t')[0] in subA_list:
		f4.write(line3)
	elif line3.strip().split('\t')[0] in subB_list:
		f5.write(line3)
	else:
		continue
f1.close()
f2.close()
f3.close()
f4.close()
f5.close()

#genome.fa
os.system("seqtk subseq %s %s > %s" % (args.sp4genome,path4 + "/subA.scaffold.name",path5 + "/subA.genome.fa"))
os.system("seqtk subseq %s %s > %s" % (args.sp4genome,path4 + "/subB.scaffold.name",path5 + "/subB.genome.fa"))

#pep
os.system("seqtk subseq %s %s> %s" % (args.sp4pep,path4 + "/subA.gene.name",path5 + "/subA.pep"))
os.system("seqtk subseq %s %s> %s" % (args.sp4pep,path4 + "/subB.gene.name",path5 + "/subB.pep"))

#cds
os.system("seqtk subseq %s %s> %s" % (args.sp4cds,path4 + "/subA.gene.name",path5 + "/subA.cds"))
os.system("seqtk subseq %s %s> %s" % (args.sp4cds,path4 + "/subB.gene.name",path5 + "/subB.cds"))

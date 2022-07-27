#import modules
import os
import argparse
import pandas as pd

#2022/07/06 created by yuanw

#print help information
parser = argparse.ArgumentParser(description='Select replicated gene pairs and construct clusters.')
parser.add_argument("-kt",help="the ks threshold of filtering ancient duplications")
parser.add_argument("-sp4gff",help="input allotetraploid species gff file")
args = parser.parse_args()

#dir
path = os.getcwd()
path1 = path + "/1.Collinear"
os.mkdir(path + '/2.Cluster')
path2 = path + "/2.Cluster"

#filter ancient duplication
os.system("cut -f 1,4 %s | awk '$2> %s {print $1}' | sed -e '1d' -e 's/|/\\n/g'|sort|uniq > %s/ancient_sp4.ks" % (path1+"/sp4_sp4.ks.txt",args.kt,path1))

f1 = open(path1 + "/ancient_sp4.ks",'r')
f2 = open(path1 + "/sp4_sp2.collinearity",'r')
f3 = open(path2 + "/sp4_sp2.collinearity.filter",'w')

ancient_list = list()
for line1 in f1.readlines():
	ancient_list.append(line1.strip())

for line2 in f2.readlines():
	if line2.startswith('#'):
		continue
	elif line2.strip().split('\t')[2] not in ancient_list:
		f3.write(line2.strip() + '\n')
	else:
		continue

f1.close()
f2.close()
f3.close()

#cluster
os.system("cut -f 1 %s/sp4_sp2.collinearity.filter|sed  -e 's/-\\s*/\\t/g' -e 's/://g' > %s/block_number.txt" % (path2,path2))
os.system("paste %s/block_number.txt %s/sp4_sp2.collinearity.filter|cut -f 1,2,4,5 > %s/sp4_sp2.collinearity.filter.num" % (path2,path2,path2))

#filter 1v2 genes
df  = pd.read_csv(path2 + "/sp4_sp2.collinearity.filter.num", sep = '\t', header = None, names=['blocknum1','blocknum2','sp2','sp4'])
dfnum = df.sp2.value_counts() #Series
one_two_list = dfnum[dfnum == 2].index.values.tolist()

f1 = open(path2 + "/sp4_sp2.collinearity.filter.num",'r')
f2 = open(path2 + "/sp4_sp2.collinearity.filter.num.1v2",'w')
for line in f1.readlines():
	if line.strip().split('\t')[2] in one_two_list:
		f2.write(line)
	else:
		continue
f1.close()
f2.close()

#format & sort (blocknum1,blocknum2,sp2,sp41,blocknum1,blocknum2,sp2,sp42)
os.system("sort -k 3 %s | sed -n '1~2p' > %s" % (path2 + "/sp4_sp2.collinearity.filter.num.1v2",path2 + "/sp4_sp2.collinearity.filter.num.1v2.sp41"))
os.system("sort -k 3 %s | sed -n '2~2p' > %s" % (path2 + "/sp4_sp2.collinearity.filter.num.1v2",path2 + "/sp4_sp2.collinearity.filter.num.1v2.sp42"))
os.system("paste %s %s |sort -n -k 1 -k 2 > %s" % (path2 + "/sp4_sp2.collinearity.filter.num.1v2.sp41",path2 + "/sp4_sp2.collinearity.filter.num.1v2.sp42",path2 + "/sp4_sp2.collinearity.filter.num.1v2.all.sort"))

#filter sp4_1 and sp4_2 in the same scaffold & blocks which <5 pairs
df2  = pd.read_csv(path2 + "/sp4_sp2.collinearity.filter.num.1v2.all.sort", sep = '\t', header = None, names=['blocknum1_1','blocknum2_1','sp2_1','sp4_1','blocknum1_2','blocknum2_2','sp2_2','sp4_2'])
dfnum2 = df2.blocknum1_1.value_counts() #Series
more5gene_list = dfnum2[dfnum2 > 5].index.values.tolist()

f1 = open(args.sp4gff,"r")
f2 = open(path2 + "/sp4_sp2.collinearity.filter.num.1v2.all.sort","r")
f3 = open(path2 + "/sp4_sp2.collinearity.filter.num.1v2.all.sort.final","w")

sp4_gff_dict = {}
for line1 in f1.readlines():
	k = line1.strip().split('\t')[1]
	v = line1.strip().split('\t')[0]
	sp4_gff_dict[k] = v

for line2 in f2.readlines():
    id1 = line2.strip().split('\t')[3]
    id2 = line2.strip().split('\t')[7]
    if sp4_gff_dict[id1] != sp4_gff_dict[id2] and int(line2.strip().split('\t')[0]) in more5gene_list:
        f3.write(line2)
    else:
        continue

f1.close()
f2.close()
f3.close()

#get cluster
os.system("cut -f 1 %s|uniq|cat -n|sed 's/ //g' > %s" % (path2 + "/sp4_sp2.collinearity.filter.num.1v2.all.sort.final",path2 + "/cluster.ifo"))

f1 = open(path2 + "/cluster.ifo","r")
f2 = open(path2 + "/sp4_sp2.collinearity.filter.num.1v2.all.sort.final","r")
f3 = open(path2 + "/sp4_sp2.1v2.cluster.number","w")

cluster_dict = {}
for line1 in f1.readlines():
	k = line1.strip().split('\t')[1]
	v = line1.strip().split('\t')[0]
	cluster_dict[k]=v

for line2 in f2.readlines():
	line2 = line2.strip().split('\t')
	f3.write("cluster" + str(cluster_dict[line2[0]])+ "\t" + str(line2[2]) + "\t" + str(line2[3]) + "\t" + str(line2[7]) + "\n")

f1.close()
f2.close()
f3.close()




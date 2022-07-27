#import modules
import os
import argparse
import pandas as pd
from Bio import SeqIO
import re

#2022/07/06 created by yuanw

#print help information
parser = argparse.ArgumentParser(description='Construct phylogenetic tree within clusters.')
parser.add_argument("-sp4gff",help="input allotetraploid species gff file")
parser.add_argument("-sp4pep",help="input allotetraploid species pep file")
parser.add_argument("-sp2pep",help="input polidy species pep file")
parser.add_argument("-outpep",help="input outgroup species pep file")
parser.add_argument("-sp4cds",help="input allotetraploid species cds file")
parser.add_argument("-sp2cds",help="input polidy species cds file")
parser.add_argument("-outcds",help="input outgroup species cds file")
parser.add_argument("-thread",help="thread number")
args = parser.parse_args()

#dir
path = os.getcwd()
path2 = path + "/2.Cluster"
os.mkdir(path + '/3.Tree')
path3 = path + "/3.Tree"

#Subcluster
f1 = open(path2 + "/sp4_sp2.1v2.cluster.number","r")

cluster_list = []
for line in f1.readlines():
	cluster_list.append(line.strip().split('\t')[0])
cluster_list = list(set(cluster_list))

f1.close()

#create cluster dir
for name in cluster_list:
       dir_name = path3 + "/" + str(name) + ".tree"
       os.mkdir(dir_name)

##divide cluster
df  = pd.read_csv(path2 + "/sp4_sp2.1v2.cluster.number", sep = '\t', header = None, names=['cluster','sp2','sp4_1','sp4_2'])
group_obj = df.groupby('cluster')
for key,value in group_obj:
	value.to_csv(path_or_buf = path3 + "/" + str(key) + ".tree/" + str(key), sep = '\t',header = False,index = False)

#sp2pep
for i in cluster_list:
	os.system("cut -f 2 %s > %s" % (path3 + "/" + str(i) + ".tree/" + str(i), path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2name"))
	os.system("seqtk subseq %s %s > %s" % (args.sp2pep,path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2name",path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2pep"))

#sp2blastout
os.system("makeblastdb -in %s -dbtype prot -parse_seqids -out %s" % (args.outpep,path3 + "/outdb"))
for i in cluster_list:
	os.system("blastp -query %s -db %s -out %s -evalue 1e-10 -num_threads %s -outfmt 6 -num_alignments 5" % (path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2pep",path3 + "/outdb",path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2blast",args.thread))

#filter best blastp out
for i in cluster_list:
	df2  = pd.read_csv(path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2blast", sep = '\t', header = None, names=['sp2','out','ident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])
	df2 = df2.sort_values("bitscore", ascending=False).drop_duplicates("sp2", keep='first').reset_index(drop=True)
	df2.to_csv(path_or_buf = path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2blast.best", sep = '\t',header = False,index = False)

#get sp2,sp41,sp42,out name
for i in cluster_list:
	f1 = open(path3 + "/" + str(i) + ".tree/" + str(i) + ".sp2blast.best","r")
	f2 = open(path3 + "/" + str(i) + ".tree/" + str(i),"r")
	f3 = open(path3 + "/" + str(i) + ".tree/" + str(i) + ".allname","w")	
	outdict = {}
	for line1 in f1.readlines():
		k = line1.strip().split('\t')[0]
		v = line1.strip().split('\t')[1]
		outdict[k] = v
	for line2 in f2.readlines():
		line2 = line2.strip().split('\t')
		if line2[1] in outdict:
			out = outdict[line2[1]]
			f3.write(str(line2[1]) + "\t" + str(out) + "\t" + str(line2[2]) + "\t" + str(line2[3]) + "\n")
		else:
			continue
	
	f1.close()
	f2.close()
	f3.close()

#split file
for i in cluster_list:
	os.system("sed 's/\\t/\\n/g' %s > %s" % (path3 + "/" + str(i) + ".tree/" + str(i) + ".allname",path3 + "/" + str(i) + ".tree/" + str(i) + ".allname.1"))
	os.system("split -l 4 %s -d %sname_" % (path3 + "/" + str(i) + ".tree/" + str(i) + ".allname.1",path3 + "/" + str(i) + ".tree/"))

#id.txt
for i in cluster_list:
	os.system("ls %s >> id.txt" % (path3 + "/" + str(i) + ".tree/" + "name*"))

#match pep
os.system("cat %s %s %s > %s" % (args.sp2pep,args.sp4pep,args.outpep,path3 + "/all.pep"))
os.system("for i in `cat id.txt`; do echo \"seqtk subseq %s $i > ${i}.pep\">>matchpep.sh; done" % (path3 + "/all.pep"))
os.system("cat matchpep.sh| xargs -i -P %s bash -c {}" % (args.thread))

#pep align
os.system("for i in `cat id.txt`; do echo \"mafft --quiet --auto ${i}.pep >${i}.afa\">>mafft.sh;done")
os.system("cat mafft.sh| xargs -i -P %s bash -c {}" % (args.thread))

#pep2cds
os.system("cat %s %s %s > %s" % (args.sp2cds,args.sp4cds,args.outcds,path3 + "/all.cds"))
all_cds_dict = {}
for record in SeqIO.parse(path3 + "/all.cds", 'fasta'):
	all_cds_dict[record.id] = record.seq

def pep2cds(afafile,cdsfile):
	f1 = open(afafile,'r')
	f2 = open(cdsfile,'w')
	for line in f1.readlines():
		line = line.strip()
		if line.startswith('>'):
			pep_name = line[1:]
			f2.write('>' + str(pep_name) + '\n')
			hit = 0
		else:
			seq = all_cds_dict[pep_name]
			for site in range(len(line)):
				if line[site] == "-":
					f2.write("---")
				else:
					f2.write(str(seq[3*hit:3*hit+3]))
					hit += 1
			f2.write("\n")
	f1.close()
	f2.close()

f1 = open("id.txt","r")
for line in f1.readlines():
	a = str(line.strip()) + ".afa"
	b = str(line.strip()) + ".cds.aln"
	pep2cds(a,b)
f1.close()

#gblocks
os.system("for i in `cat id.txt`; do Gblocks ${i}.cds.aln -t=c -b4=5 -b5=a;done")
os.system("for i in `cat id.txt`; do echo \"seqkit seq -w 0 ${i}.cds.aln-gb > ${i}.cds.aln-gb.1\">>gblock.1.sh;done")
os.system("cat gblock.1.sh| xargs -i -P %s bash -c {}" % (args.thread))

#paste in cluster & name
f1 = open("cluster.dir.txt","w")
for i in cluster_list:
	f1.write(str(path3) + "/" + str(i) + ".tree/" + "\n")
f1.close()

os.system("for i in `cat cluster.dir.txt`;do echo \"paste ${i}*gb.1 > ${i}pretree.all\">>paste.sh;done")
os.system("cat paste.sh| xargs -i -P %s bash -c {}" % (args.thread))

os.system("for i in `cat cluster.dir.txt`;do sed -e 's/ //g' -e 's/\t//g' -e '1c >sp2' -e '3c >sp41' -e '5c >sp42' -e '7c >out' ${i}pretree.all > ${i}pretree.all.name;done")

f1 = open("cluster.dir.txt","r")
for line in f1.readlines():
	line = line.strip()
	os.system("mv %s %s" % (str(line) + "pretree.all.name",path3 + "/" + str(line.split('/')[-2]) + ".fa"))
f1.close()

#tree construction
f1 = open("cluster.dir.txt","r")
for line in f1.readlines():
	line = line.strip()
	os.system("iqtree -s %s -m MFP -bb 1000 -T %s -o out" % (path3 + "/" + str(line.split('/')[-2]) + ".fa",args.thread))
f1.close()

#statistic
#summary
f1 = open(path3 + "/all.treefile.bootstrap","w")
for i in cluster_list:
	f2 = open(path3 + "/" + i + ".tree.fa.treefile","r")
	s = f2.read().strip()
	if s.startswith("((sp2"):
		number = re.findall('([\d+\.]+)',s)
		for item in number:
			s = s.replace(':' + item, '')
		s = s.replace(number[4], '')
		f1.write(str(i) + "\t" + str(s) +  "\t" + str(number[4])  + "\n")
	if s.startswith("(sp2"):
		number = re.findall('([\d+\.]+)',s)
		for item in number:
			s = s.replace(':' + item, '')
		s = s.replace(number[6], '')
		f1.write(str(i) + "\t" + str(s) +  "\t" + str(number[6])  + "\n")
	f2.close()
f1.close()

#remove
os.system("rm -rf *.txt")
os.system("rm -rf *.sh")

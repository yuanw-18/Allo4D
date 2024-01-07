#import modules
import os
import argparse
from Bio import SeqIO
import pandas as pd

#2022/07/06 created by yuanw

#print help information
parser = argparse.ArgumentParser(description='Align allotetraploid and diploid genome data and obtain their collinear relationship.')
parser.add_argument("-sp4pep",help="input allotetraploid species pep file")
parser.add_argument("-sp4cds",help="input allotetraploid species cds file")
parser.add_argument("-sp4bed",help="input allotetraploid species bed file")
parser.add_argument("-sp2pep",help="input diploid species pep file")
parser.add_argument("-sp2bed",help="input diploid species bed file")
parser.add_argument("-thread",help="thread number")
args = parser.parse_args()

##step1:Align allotetraploid and diploid genome data and obtain their collinear relationship.
path = os.getcwd()
os.mkdir(path + '/1.Collinear')
path1 = path + "/1.Collinear"

#blastp
os.system("makeblastdb -in %s -dbtype prot -parse_seqids -out %s/%s.db" % (args.sp4pep,path1,args.sp4pep))
#sp4&sp2 blastp
os.system("blastp -query %s -db %s/%s.db -out %s/sp4_sp2.blast -evalue 1e-10 -num_threads %s -outfmt 6 -num_alignments 5" % (args.sp2pep,path1,args.sp4pep,path1,args.thread))
#sp4&sp4 blastp
os.system("blastp -query %s -db %s/%s.db -out %s/sp4_sp4.blast -evalue 1e-10 -num_threads %s -outfmt 6 -num_alignments 5" % (args.sp4pep,path1,args.sp4pep,path1,args.thread))

#filter blastp result(identity>30%,coverage>30%)
os.system("cat %s %s > sp4_sp2.pep" % (args.sp4pep,args.sp2pep))
peplen_dict = {}

for record in SeqIO.parse("sp4_sp2.pep", 'fasta'):
	peplen_dict[record.id] = len(record.seq)

def FilterBlastp(blastout,filterout):
	f1 = open(blastout,'r')
	f2 = open(filterout,'w')
	for line in f1.readlines():
		line = line.strip().split('\t')
		ident = float(line[2])
		aln_len = line[3]
		len1 = peplen_dict[str(line[0])]; len2 = peplen_dict[str(line[1])]
		cov1 = int(aln_len)/int(len1); cov2 = int(aln_len)/int(len2)
		if ident > 30 and cov1 > 0.3 and cov2 > 0.3:
			f2.write("\t".join(line) + '\n')
		else:
			continue
	f1.close()
	f2.close()
	return

FilterBlastp(path1 + "/sp4_sp2.blast",path1 + "/sp4_sp2.blast.filter")
FilterBlastp(path1 + "/sp4_sp4.blast",path1 + "/sp4_sp4.blast.filter")

#MCscan

#blastp prepare
sp4_sp2bfilter = path1 + "/sp4_sp2.blast.filter"
sp4_sp2bfinal = path1 + "/sp4_sp2.blast"

sp4_sp4bfilter = path1 + "/sp4_sp4.blast.filter"
sp4_sp4bfinal = path1 + "/sp4_sp4.blast"

os.system("mv %s %s" % (sp4_sp2bfilter,sp4_sp2bfinal))
os.system("mv %s %s" % (sp4_sp4bfilter,sp4_sp4bfinal))

#bed prepare
os.system("cat %s %s > %s/sp4_sp2.bed" % (args.sp4bed,args.sp2bed,path1))
os.system("cat %s %s > %s/sp4_sp2.gff" % (args.sp4bed,args.sp2bed,path1))
os.system("cp %s %s/sp4_sp4.bed" % (args.sp4bed,path1))
os.system("cp %s %s/sp4_sp4.gff" % (args.sp4bed,path1))

#run MCscan
os.system("MCScanX %s/sp4_sp2" % (path1))
os.system("MCScanX %s/sp4_sp4" % (path1))

#filter ancient duplication(ks within sp4)
os.system("grep -v '#' %s/sp4_sp4.collinearity|cut -f 2,3|sort|uniq > %s/sp4_sp4.collinearity.homolog" % (path1,path1))

#ks calculation
#pep,cds (filter cds:pep != 3)
sp4_cds_dict = {}
for record in SeqIO.parse(args.sp4cds, 'fasta'):
	sp4_cds_dict[record.id] = record.seq

sp4_pep_dict = {}
for record in SeqIO.parse(args.sp4pep, 'fasta'):
	sp4_pep_dict[record.id] = record.seq

f1 = open(path1 + "/sp4_sp4.collinearity.homolog",'r')
f2 = open(path1 + "/sp4_sp4.collinearity.homolog.filter",'w')
for line in f1.readlines():
	name1 = line.strip().split("\t")[0]
	name2 = line.strip().split("\t")[1]
	rate1 = len(sp4_cds_dict[name1])/len(sp4_pep_dict[name1])
	rate2 = len(sp4_cds_dict[name2])/len(sp4_pep_dict[name2])
	if rate1==3 and rate2==3:
		f2.write(str(name1) + "\n" + str(name2) + "\n")
	else:
		continue
f1.close()
f2.close()
#pep match (sp4_pep_dict created before)
f1 = open(args.sp4pep,'r')
f2 = open(path1 + "/sp4_sp4.collinearity.homolog.filter",'r')
f3 = open(path1 + "/sp4_sp4.collinearity.homolog.filter.pep",'w')
for line in f2.readlines():
	line = line.strip()
	f3.write(">" + str(line) + "\n" + str(sp4_pep_dict[line]) + "\n")
f1.close()
f2.close()
f3.close()
#pep align
homolog_filter_pep = path1 + "/sp4_sp4.collinearity.homolog.filter.pep"
pepsh = "cat " + homolog_filter_pep + " |xargs -n 4 -P " + args.thread + " bash -c " + "\'" + 'echo -e \"$0\\n$1\\n$2\\n$3\" | mafft --quiet -' + "\'" + ">" + path1 + "/all.pep.aln"
f1 = open(path1 + "/mafft.sh",'w')
f1.write(pepsh)
f1.close()
os.system("sh %s/mafft.sh" % (path1))
#pep2cds
f1 = open(path1 + "/all.pep.aln",'r')
f2 = open(path1 + "/all.cds.aln",'w')
for line in f1.readlines():
    line = line.strip()
    if line.startswith('>'):
        pep_name = line[1:]
        f2.write('>' + str(pep_name) + '\n')
        hit = 0
    else:
        seq = sp4_cds_dict[pep_name]
        for site in range(len(line)):
            if line[site] == "-":
                f2.write("---")
            else:
                f2.write(str(seq[3*hit:3*hit+3]))
                hit += 1
        f2.write("\n")
f1.close()
f2.close()
#afa2axt
os.system("seqkit seq -w 0 %s > %s" % (path1 + "/all.cds.aln",path1 + "/all.cds.aln.1"))

f1 = open(path1 + "/all.cds.aln.1",'r')
f2 = open(path1 + "/all.cds.aln.axt",'w')
id_lst = []
seq_lst = []

for line in f1.readlines():
	line = line.strip()
	if line.startswith('>'):
		line = line[1:]
		id_lst.append(line)
	else:
		seq_lst.append(line)

for i in range(0,len(id_lst),2):
	j = i+1
	f2.write(str(id_lst[i]) + '|' + str(id_lst[j]) + '\n' +
			str(seq_lst[i]) + '\n' + str(seq_lst[j]) + '\n' + '\n')

f1.close()
f2.close()
#ks
os.system("KaKs_Calculator -i %s -o %s -m YN" % (path1 + "/all.cds.aln.axt",path1 + "/sp4_sp4.ks.txt"))
#ks.density.fig
df = pd.read_csv(path1 + "/sp4_sp4.ks.txt",header=0,sep='\t')
density = df.Ks.dropna().plot.hist(density=True,bins=1000,xlim=(0,5),xticks=(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5))
fig = density.get_figure()
fig.savefig(path1 +'Ks_density.png')

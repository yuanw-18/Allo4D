# Allo4D (Allotetraploid Dividing)
Allo4D is a pipeline to distinguish allotetraploid subgenomes based scaffold level. 

# Dependencies
**1.blastp**
**2.MCScanX**
**3.mafft**
**4.KaKs_Calculator**
**5.iqtree**
**6.seqkit**
**7.seqtk**
**8.biopython**
**9.pandas**
**10.Gblocks**
```
conda install -c bioconda blast mafft kakscalculator2 iqtree seqkit seqtk
conda install -c conda-forge biopython
conda install -c anaconda pandas
conda install bioconda::gblocks
```
MCScanX: https://github.com/wyp1125/MCScanX

# Input data
**allotetraploid**：sp4.pep;sp4.cds;sp4.bed;sp4.genome.fa

**relative diploidy**：sp2.pep;sp2.cds;sp2.bed

**outgroup**：outgroup.pep;outgroup.cds


For example data, we obtained four genomes (_Oryza sativa_ (_indica_ group), _Oryza sativa_ (_japonica_ group), _Oryza punctata_, _Leersia perrieri_) in RiceRelativesGD (http://ibi.zju.edu.cn/ricerelativesgd/), and genomes were cut according to the criterion of N50 = 1Mb. The genomes of _Oryza sativa_ (_indica_ group) and _Oryza punctata_ were combined as a putative allotetraploid, _Oryza sativa_ (_japonica_ group) as relative diploid, and _Oryza punctata_ as outgroup.

Be careful, bed file is only four columns: 
```
#scaffold gene start end
R498Chr1.ctg1   indica_OsR498G0100000700.01.T01 17040   23888
R498Chr1.ctg1   indica_OsR498G0100001000.01.T01 24847   26162
R498Chr1.ctg1   indica_OsR498G0613293000.01.T01 28450   28624
R498Chr1.ctg1   indica_OsR498G0100001100.01.T01 31518   32476
R498Chr1.ctg1   indica_OsR498G0100001500.01.T01 36552   39342
R498Chr1.ctg1   indica_OsR498G0100001700.01.T01 40581   43313
R498Chr1.ctg1   indica_OsR498G0100002000.01.T01 44826   46136
R498Chr1.ctg1   indica_OsR498G0100002100.01.T01 47224   50422
R498Chr1.ctg1   indica_OsR498G0100002300.01.T01 51254   52406
R498Chr1.ctg1   indica_OsR498G0100002500.01.T01 52877   53327
```

genome file, pep file and cds file are common fasta format.

# How to run Allo4D ?

## Step1: Align allotetraploid and diploid genome data and obtain their collinear relationship.
```
python Allo4Dalign.py -thread 40 \
 -sp4pep sp4.pep -sp2pep sp2.pep \
 -sp4bed sp4.bed -sp2bed sp2.bed \
 -sp4cds sp4.cds
```
### Output files (You can check in ./1.Collinear/)
```
  75M Jul  8 14:27 all.cds.aln
  75M Jul  8 14:34 all.cds.aln.1
  75M Jul 10 09:37 all.cds.aln.axt
  27M Jul  7 22:20 all.pep.aln
 290K Jul 26 15:32 ancient_sp4.ks
  210 Jul  7 22:11 mafft.sh
 5.1M Jul  7 14:49 sp4.pep.db.phr
 623K Jul  7 14:49 sp4.pep.db.pin
 312K Jul  7 14:49 sp4.pep.db.pog
 5.6M Jul  7 14:49 sp4.pep.db.psd
 108K Jul  7 14:49 sp4.pep.db.psi
  29M Jul  7 14:49 sp4.pep.db.psq
  12M Jul  7 16:19 sp4_sp2.blast
 3.3M Jul  7 16:57 sp4_sp2.collinearity
 6.7M Jul  7 16:47 sp4_sp2.bed
  60K Jul  7 16:57 sp4_sp2.html
  28M Jul  7 16:19 sp4_sp4.blast
 2.0M Jul  7 16:57 sp4_sp4.collinearity
 1.4M Jul  7 17:21 sp4_sp4.collinearity.homolog
 1.4M Jul  7 18:50 sp4_sp4.collinearity.homolog.filter
  24M Jul  7 21:26 sp4_sp4.collinearity.homolog.filter.pep
 4.6M Jul  7 16:50 sp4_sp4.bed
  36K Jul  7 16:57 sp4_sp4.html
 5.6M Jul 10 09:42 sp4_sp4.ks.txt
 568K Jul  7 16:57 sp4_sp4.tandem
 8.4K Jul 22 18:16 Ks_density.png
```

## Step2: Select replicated gene pairs and construct clusters.
Find the ancient Ks threshold according to the 1.Collinear/Ks_density.png. 
For example data,Ks threshold = 0.5
```
python Allo4Dcluster.py -kt 0.5 -sp4bed sp4.bed
```
### Output files (You can check in ./2.Cluster/)
```
212K Jul 26 15:32 block_number.txt
4.0K Jul 26 15:32 cluster.ifo
977K Jul 26 15:32 sp4_sp2.1v2.cluster.number
2.1M Jul 26 15:32 sp4_sp2.collinearity.filter
1.8M Jul 26 15:32 sp4_sp2.collinearity.filter.num
1.4M Jul 26 15:32 sp4_sp2.collinearity.filter.num.1v2
1.4M Jul 26 15:32 sp4_sp2.collinearity.filter.num.1v2.all.sort
1.3M Jul 26 15:32 sp4_sp2.collinearity.filter.num.1v2.all.sort.final
706K Jul 26 15:32 sp4_sp2.collinearity.filter.num.1v2.sp41
629K Jul 26 15:32 sp4_sp2.collinearity.filter.num.1v2.sp42
```
sp4_sp2.collinearity.filter.num file: Please make sure the diploid gene in the third column, and the allotetraploid gene in the fourth column. If the order is reversed, amend line 45 in Allo4Dcluster.py.
```
os.system("paste %s/block_number.txt %s/sp4_sp2.collinearity.filter|awk -v FS='\t' -v OFS='\t' '{print $1,$2,$5,$4}' > %s/sp4_sp2.collinearity.filter.num" % (path2,path2,path2))
```


## Step3: Construct phylogenetic tree within clusters.
```
python Allo4Dtree.py -thread 40 \
 -sp2pep sp2.pep -sp4pep sp4.pep -outpep out.pep \
 -sp2cds sp2.cds -sp4cds sp4.cds -outcds out.cds
```
### Output files (You can check in ./3.Tree/)
Cluster and tree files.

## Setp4:Divide subgenomes.
```
python Allo4Ddivide.py \
 -sp4bed sp4.bed -sp4genome sp4.genome.fa \
 -sp4pep sp4.pep -sp4cds sp4.cds
```
### Output files (You can check in ./4.groupAB/ and ./5.result/)
### ./4.groupAB/
```
 403K Jul 26 17:25 groupAB.gene.name
  17K Jul 26 17:25 groupAB.gene.name.count
 1.3M Jul 26 17:25 groupAB.result
  19K Jul 26 17:25 groupAB.summary
 5.1K Jul 26 17:25 group_A_B.txt
   33 Jul 26 17:25 group_B_A.txt
   77 Jul 26 17:25 group_NS.txt
  19K Jul 26 17:25 group.txt
  13K Jul 26 17:25 scaffold.name
  20K Jul 26 17:25 sp4.genome.fa.len
 1.1M Jul 26 17:25 subA.gene.name
 5.1K Jul 26 17:25 subA.scaffold.name
 972K Jul 26 17:25 subB.gene.name
 7.3K Jul 26 17:25 subB.scaffold.name
```

### ./5.result/ (the final result)
```
  35144893 Jul 26 17:25 subA.cds
 343391956 Jul 26 17:25 subA.genome.fa
   2157791 Jul 26 17:25 subA.bed
  12521059 Jul 26 17:25 subA.pep
  52915871 Jul 26 17:25 subB.cds
 370110443 Jul 26 17:25 subB.genome.fa
   2329783 Jul 26 17:25 subB.bed
  18354505 Jul 26 17:25 subB.pep
```
# Contact me
Wang Yuan: wangy959@mail3.sysu.edu.cn (E-mail can be in Chinese)

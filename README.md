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
```
conda install -c bioconda blast mafft kakscalculator2 iqtree seqkit seqtk
conda install -c conda-forge biopython
conda install -c anaconda pandas
```
MCScanX: https://github.com/wyp1125/MCScanX

# Input data
**allotetraploid**：sp4.pep;sp4.cds;sp4.gff;sp4.genome.fa

**relative diploidy**：sp2.pep;sp2.cds;sp2.gff

**outgroup**：outgroup.pep;outgroup.cds


For example data, we obtained four genomes (_Oryza sativa_ (_indica_ group), _Oryza sativa_ (_japonica_ group), _Oryza punctata_, _Leersia perrieri_) in RiceRelativesGD (http://ibi.zju.edu.cn/ricerelativesgd/), and genomes were cut according to the criterion of N50 = 1Mb. The genomes of _Oryza sativa_ (_indica_ group) and _Oryza punctata_ were combined as a putative allotetraploid, _Oryza sativa_ (_japonica_ group) as relative diploid, and _Oryza punctata_ as outgroup.

Be careful, gff file is only four columns: scaffold name; gene name; start position; end position
```
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
# How to run Allo4D ?

## Step1: Align allotetraploid and diploid genome data and obtain their collinear relationship.
```
python Allo4Dalign.py -thread 40 \
 -sp4pep sp4.pep -sp2pep sp2.pep \
 -sp4gff sp4.gff -sp2gff sp2.gff \
 -sp4cds sp4.cds
```

## Step2: Select replicated gene pairs and construct clusters.
```
#find the ancient ks threshold according to the 1.Collinear/Ks_density.png. For example data,ks threshold = 0.5
python Allo4Dcluster.py -kt 0.5 -sp4gff sp4.gff
```

## Step3: Construct phylogenetic tree within clusters.
```
python Allo4Dtree.py -thread 40 \
 -sp2pep sp2.pep -sp4pep sp4.pep -outpep out.pep \
 -sp2cds sp2.cds -sp4cds sp4.cds -outcds out.cds
```

## Setp4:Divide subgenomes.
```
python Allo4Ddivide.py \
 -sp4gff sp4.gff -sp4genome sp4.genome.fa \
 -sp4pep sp4.pep -sp4cds sp4.cds
```

# Output files (You can check in 5.result/)
```
  35144893 Jul 26 17:25 subA.cds
 343391956 Jul 26 17:25 subA.genome.fa
   2157791 Jul 26 17:25 subA.gff
  12521059 Jul 26 17:25 subA.pep
  52915871 Jul 26 17:25 subB.cds
 370110443 Jul 26 17:25 subB.genome.fa
   2329783 Jul 26 17:25 subB.gff
  18354505 Jul 26 17:25 subB.pep
```
# Contact me
Wang Yuan: wangy959@mail2.sysu.edu.cn (E-mail can be in Chinese)

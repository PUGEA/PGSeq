# Annotation Type

* * *

PGSeq supports three annotation databases: RefSeq, Ensembl and UCSC. All annotation information from three databases are downloaded from the [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables?command=start). And the annotation information from Ensembl has  an [alternative (Ensembl website)](http://asia.ensembl.org/index.html). Therefore, there are   four types of annotation files.

The four types are:

1.  refGene2.  ensGene3.  knownGene4.  Ensembl

The first three types of annotation files are downloaded from UCSC Table Browser. The 'Ensembl' type of annotation file is downloaded from Ensembl website. 

The 'ensGene' and 'Ensembl' have the same annotation information as Ensembl database, but they have the differential format.

&nbsp;

* * *

### In the follow section, we will give four examples for four types of the annotation file respectively.

&nbsp;

### 1. refGene

*   refGene.txt is downloaded from the UCSC Table Browser. (genome: Human, track: RefSeq Genes, table: refGene)
<table width="100%" border="0">
  <tr>
    <th align="left" bgcolor="#CCCCCC" scope="col">

$ python ./PGSeq/preprocessAlignment<span class="prettyprint">.py</span> -t 1 -a refGene.txt -o HumPrefix -d data1.sam,data2.sam,data3.sam
</th>
  </tr>
</table>

### &nbsp;

### 2. ensGene

*   ensGene.txt is downloaded from  the  UCSC Table Browser. (genome: Human, track: Ensembl Genes, table: ensGene)
<table width="100%" border="0">
  <tr>
    <th align="left" bgcolor="#CCCCCC" scope="col">

$ python ./PGSeq/preprocessAlignment<span class="prettyprint">.py</span> -t 2 -a ensGene.txt -o HumPrefix -d data1.sam,data2.sam,data3.sam
</th>
  </tr>
</table>

### &nbsp;

### 3. knownGene

*   In order to obtain the relationship between genes and isoforms, PGSeq need two annotation files: knownGene.txt and knownIsoform.txt*   knownGene.txt is downloaded from   the UCSC Table Browser. (genome: Human, track: UCSC Genes, table: knownGene)
*   knownIsoform.txt is downloaded from the UCSC Table Browser. (genome: Human, track: UCSC Genes, table: knwonIsoforms)
<table width="100%" border="0">
  <tr>
    <th align="left" bgcolor="#CCCCCC" scope="col">

$ python ./PGSeq/preprocessAlignment<span class="prettyprint">.py</span> -t 3 -a knownGene.txt,knownIsoforms.txt -o HumPrefix -d data1.sam,data2.sam,data3.sam
</th>
  </tr>
</table>

&nbsp;

### 4. Ensembl

*   Homo_sapiens.GRCH37.67.gtf is downloaded from the Ensembl website. Because we only download the 'gtf' format, this type is diffenent from the above three tpyes.
<table width="100%" border="0">
  <tr>
    <th align="left" bgcolor="#CCCCCC" scope="col">

$ python ./PGSeq/preprocessAlignment<span class="prettyprint">.py</span> -t 4 -a Homo_sapiens.GRCH37.67.gtf -o HumPrefix -d data1.sam,data2.sam,data3.sam
</th>
  </tr>
</table>

Since the 'ensGene' and 'Ensembl' have the same annotation information from Ensembl database, we suggest that users select the 'ensGene' type.

&nbsp;

</body>
</html>
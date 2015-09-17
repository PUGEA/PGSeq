# PGSeq
A bioinformatics tool for estimating gene and isoform expression level


Xuejun Liu (xuejun.liu@nuaa.edu.cn)


* * *

## <a name="introduction"></a> Introduction

PGSeq   uses the Poisson model to fit read counts for each gene and 
  use gene-specific Gamma-distributed latent variables to capture the variability of read 
  sequencing preference for every exonic position. The bias property modeled in our 
  model is shared across all conditions for each individual gene, and automatically 
  captures all intrinsic position-specific effects.

PGSeq pacakage is the software for estimating gene and isoform expression levels from RNA-seq data, and provides a user-friendly interface and supports parallel computation.

* * *

## <a name="compilation"></a>  Installation

You can click [here](https://github.com/PUGEA/PGSeq/tree/master/PGSeq.1.1) to download the PGSeq software. 

To compile PGSeq in Linux, run in the PGSeq folder.

```shell
$ bash  setup.sh
```

### Requirements:

*   PGSeq  implementation uses Python (v.2.7) to pre-process the RNA-seq data and C language to calculate the gene and isoform expression levels.
*   In PGSeq, the Python codes use two special modules, [NumPy](http://www.numpy.org/) and [PP](http://www.parallelpython.com/) (parallel python).
*   PGSeq uses Bowtie2 to align reads to transcript reference  sequences, so you need to have [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) installed.

&nbsp;

* * *

## <a name="usage"></a> Usage

### Step 1. Aligning Reads

Create index and align reads to reference transcript sequences:

```shell
$ bowtie2 -t -p 4 -q -k 100 -x ensGene.ref_transcript.index -1 read_1.fq -2 read_2.fq -S data1.sam --no-hd  --no-unal --no-mixed --no-discordant
```

*   PGSeq  uses the transcript reference sequences, which can be downloaded from [UCSC](http://genome.ucsc.edu/) or [Ensembl](http://asia.ensembl.org/index.html) website.
*   For paired-end reads, PGseq uses options '--no-mixed' and '--no-discordant'.
*   For the SAM format, we use options, '--no-hd' and '-no-unal', to avoid the header and reads that fail to align.
*   If you want to get more  information for the usage of Bowtie2, please click [here](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) to visit the Bowtie2 website.

### &nbsp;

### Step 2. Pre-processing Alignment Files

PGSeq needs to pre-process the alignment files using the following steps, which include pre-computing the probabilities for each alignments and extracting the count data for each gene:

```shell
$ python ./PGSeq/preprocessAlignment.py -t 2 - a ensGene.txt -o EnsPrefix -d data1.sam,data2.sam,data3.sam
```
or
```shell
$ python ./PGSeq/preprocessAlignment.py -t 2 - a ensGene.txt -o EnsPrefix -d data1.sam,data2.sam,data3.sam    -s selected_genes
```

Options:

*   -t/--AnnotationType &lt;int&gt;: an integer to select the type of annotation. refGene:1, ensGen:2, knownGene:3 and Ensmbl: 4.  For detailed interpretation, please click [here](https://github.com/PUGEA/PGSeq/blob/master/Doc/README.md).
*   -a/--AnnotationFile &lt;path&gt;: the annotation file, eg; ensGene.txt. For detailed interpretation, please click [here](https://github.com/PUGEA/PGSeq/blob/master/Doc/README.md).
*   -o/--AnnotationPrefix &lt;text&gt;: the header of  temporary files which include differential annotation information. Users can set arbitrary text, eg.: 'EnsPrefix' or 'abcd'.
*   -d/--AlignmentFiles &lt;path&gt; : input  all alignment files.  eg: data1.sam,data2.sam.
*   -s/--SelectedGenes &lt;path&gt;: optional parameters, file containing gene names. If you only calculate a subset of genes/isoforms, you can choose this option. Otherwise, PGSeq will calculate all genes/isoforms in the annotation file.

### &nbsp;

### Step 3. Calculating Expression Values

PGSeq starts calculating expression values using the following commands:

```shell
$ python ./PGSeq/calculateExpression.py -a EnsPrefix --log
```
or

```shell
$ python ./PGSeq/calculateExpression.py -a EnsPrefix -s selected_genes --log
```

Options:

*   -o/--AnnotationPrefix &lt;text&gt;: the text must be the same as the '--AnnotationPrefix' in the Step 2.
*   -s/--SelectedGenes &lt;path&gt;: optional parameters, file contianing gene names. If you only calculate a subset of genes/isoforms, you can choose this option. Otherwise, PGSeq will calculate all genes/isoforms in the annotation file.
*   -l/--log:  optional parameters. This option will calculate the logged values of gene/isoform expression.

Output Files:

The PGSeq will produce four output files, which include genes.mean, genes.std, isoforms.mean and isoforms.std.

 If you choose  '-l/--log', the PGSeq will produce the additional four output files, which include genes.mean.log, genes.std.log, isoforms.mean.log and isoforms.std.log.

Description of output files:

*   **genes/isoforms.mean:** these files contain the gene or isoform expression. The first column is the name of gene or isoform and the remaining columns are gene or isoform expression for each alignment file.
*   **genes/isoforms.std**: these files contain the standard deviation of estimated gene or isoform expression. The first column is the name of gene or isoform and the remaining columns are the standard deviation of gene or isoform expression for each alignment file.
*   **genes/isoforms.mean.log**: same as genes/isoforms.mean files, but contain logarithmic expression. 
*   **genes/isoforms.std.log**: same as genes/isoforms.std files, but contain the standard deviation of the estimated logarithmic expression.

&nbsp;

* * *

## <a name="example"></a> Example

Here, we give an example of the usage of PGSeq. 

*   Annotation file: [ensGene.txt ](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=404134387_Acyi8a1auSUTdR5AeBjiGWTPAY3w&amp;clade=mammal&amp;org=Human&amp;db=hg19&amp;hgta_group=genes&amp;hgta_track=ensGene&amp;hgta_table=0&amp;hgta_regionType=genome&amp;position=chr21%3A33%2C031%2C597-33%2C041%2C570&amp;hgta_outputType=primaryTable&amp;hgta_outFileName=)
*   Reference sequence file: [ensGene.ref_transcript.fasta](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=404134387_Acyi8a1auSUTdR5AeBjiGWTPAY3w&amp;clade=mammal&amp;org=Human&amp;db=hg19&amp;hgta_group=genes&amp;hgta_track=ensGene&amp;hgta_table=0&amp;hgta_regionType=genome&amp;position=chr21%3A33%2C031%2C597-33%2C041%2C570&amp;hgta_outputType=primaryTable&amp;hgta_outFileName=).
*   Sequenced reads: the paired-end simulation dataset from the following section.
*   Selected Genes: [selected.genes.info](https://github.com/PUGEA/PGSeq/blob/master/Example.Data/selectedgenes.zip).

You can use the following commands to run the  example.

```shell
$bowtie2 -t -p 4 -q -k 100 -x ensGene.ref_transcript.index -1 read_1.fq -2 read_2.fq -S data1.sam --no-hd  --no-unal --no-mixed --no-discordant ( The same treatment for the rest datasets.)
 
$AlignData='data1.sam,data2.sam,data3.sam,data4.sam,data5.sam,data6.sam,data7.sam';

$python ./PGSeq/preprocessAlignment.py -t 2 - a ensGene.txt -o EnsPrefix -d $AligData -s  selected.genes.info

$python ./PGSeq/calculateExpression.py -a EnsPrefix -s  selected.genes.info --log
```

&nbsp;

* * *

## <a name="example" id="example"></a> Simulation Data 

 We generated  two simulated datasets using our model based on the calculated parameters from the qRT-PCR validated genes of HBR sample in MAQC dataset. The two simulation datasets both used the Ensembl (NCBI37/hg19) as a reference.

*   [Single-end simulation dataset](https://github.com/PUGEA/PGSeq/blob/master/Example.Data/SESimulationData.zip). (35bp, 7 technical replicates, about 100 million reads for each replicate)
*   [Paired-end simulation dataset](https://github.com/PUGEA/PGSeq/blob/master/Example.Data/PESimulationData.zip) (50bp, 7 technical replicates, about 100 million reads for each replicate, fragment length: N(206, 19.6))

&nbsp;

* * *

## <a name="authors"></a> Authors

The PGSeq algorithm is developed by Xuejun Liu and Li Zhang. The PGSeq software is mainly implemented by  Li Zhang.

* * *

</body></html>

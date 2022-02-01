# DATing Evolutionary Events and Divergence times (DATED)

DATED provides an efficient single step solution to estimate the level of synonymous substitution (Ks) between paralogous and orthologous sequence pairs. The software utilizes multiprocessing library to speedup Ks calculation for input sequence pairs. The output from the pipeline can be used in the `divergence.R` script to perform mixture model analysis of Ks distribution.



### Installation

Create and activate a Conda environment named Ks

`conda create --name ks python=3.5` <br>
`conda activate ks`

Downloading and installing software: <br>
1.	DATED <br>
`git clone https://github.com/ChuShin/dated`
2.  ClustalW <br>
`conda install -c bioconda clustalw`
4.  PAL2NAL <br>
`conda install -c bioconda pal2nal`
4.	PAML <br>
`conda install -c bioconda paml`


### Identification of Paralogs:
An all-against-all protein sequence similarity (BLASTP with E-value, high-scoring segment pair (HSP) length and sequence identify cut-offs) search can be used to identify paralogous genes within a plant species for which completely annotated genome sequence is available. In the absence of a completely annotated genome sequence, transcript sequences assembled from RNA-Seq data can be used to identify homologs. In such case, the open reading frame for each transcript has to be predicted and corresponding translated amino acid sequence should be deduced.

### Identification of Orthologs:
Reciprocal best blast hit method can be used to detect orthologous genes between two related species. 


### Usage
`dated.py pep.fa cds.fa blast_pairlist > blast_pairlist.ks`


### Interpreting the results


### Citation






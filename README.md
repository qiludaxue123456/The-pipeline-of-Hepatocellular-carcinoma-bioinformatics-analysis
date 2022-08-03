# The pipeline of Hepatocellular carcinoma bioinformatics analysis 

In the study “An Adenylate Uridylate (AU)-Rich Element Genes Related Prognosis Model for Hepatocellular Carcinoma based on Bioinformatics”, the RNA-seq data and clinical information of HCC patients were acquired from The Cancer Genome Atlas (TCGA) database and Gene Expression Omnibus (GEO) database. A total of 189 differentially expressed Adenylate Uridylate (AU)-Rich Element genes (DE-AREGs) between normal and HCC samples were identified in HCC. Prognostic related DE-AREGs were selected by univariate Cox regression analysis and least absolute shrinkage and selection operator (LASSO) cox analysis, including CENPA, TXNRD1, RABIF, UGT2B15, and SERPINE1. Furthermore, multivariate Cox regression analysis was used to construct an Adenylate Uridylate (AU)-Rich Element genes (AREGs) related signature. Moreover, Kaplan-Meier curves and the receiver operating characteristic (ROC) curves revealed that the AREGs-related signature could effectively predict the prognosis of HCC patients. Additionally, T stage and risk score were screened as independent prognostic factors and C-index of 1, 3 and 5 years indicated that the nomogram performed well. Functional analysis revealed that 2242 Gene Ontology (GO) and 72 Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways were enriched in the high-risk group, but the low-risk group was enriched for only 2 GO terms and 13 KEGG pathways. Immune-related analyses indicated that T cell and B cell receptor abundance, microvascular endothelial cells (MVE), lymphatic endothelial cells (lye), pericytes, and stromal cells and the six immune checkpoints were significant differences between the high- and low-risk groups. In conclusion, an AREGs-related signature based on five DE-AREGs was constructed and could act as a prognostic indicator of HCC patients.
Here, we will provide the steps of these analyses and some of the key document information used in these analyses.

## Data downloading and processing

### TCGA data downloading

TCGA data sources, retrieve the TCGA – LIHC from GDC official website https://portal.gdc.cancer.gov/repository, and add the FPKM file of gene expression profile to cart, then download the cart files and metadata files directly, the metadata file downloaded at that time can be found in attachment Metadata.cart.2021-12-15.json, the downloaded data is then consolidated into an expression matrix using putFilestoonedir.pl and mrna_merge.pl scripts。

``` perl
# go to the directory where the expression data is decompressed
perl path/to/putFilestoonedir.pl #put all gz files in one directory
# go to the folder where the data was placed in the previous step
perl path/to/mrna_merge.pl path/to/Metadata.cart.2021-12-15.json #merge the expression data of a single sample into a matrix
```

after that, through ensemblToSymbol.pl script in combination with Homo_sapiens.GRCh38.104.chr.gtf (http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/) annotation files, the gene ID in the matrix were converted to obtain mRNA.symbol.txt finally. 

``` perl
perl path/to/ensemblToSymbol.pl path/to/Homo_sapiens.GRCh38.104.chr.gtf path/to/ensemblmatrix.txt path/to/symbolmatrix.txt #convert ID
```

Similarly, search TCGA-LIHC from GDC official website, add clinical data to CART, download the cart file directly, and then get the clinical information via getClinical.pl.

``` perl
# enter the directory where the clinical information is decompressed
perl path/to/getClinical.pl  #extract clinical information
```

### GEO data downloading

GEO data sources, download GSE14520-gPL3921_series_matrix.txt.gz data that use GPL3921 platform from GEO official website (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/), extract the expression matrix and obtain the matrix file GSE14520_probeExpr. TXT, at the same time, download the GPL3921 platform annotation file gPL3921-25447.txt, then the probe annotation file gPL3921-25447 probe-symble.TXT was obtained. 

At the same time, download the supplementary file gse14520_extra_supplements.txt. gz of this data set and extract clinical information to obtain gse14520_clinical.txt.

## The above is the initial input file of the data analysis script code.R, and all subsequent analysis is based on the above data.

Install R and RSrudio software, then open RSrudio software and execute the code in code.R by line.

## Other key Information

Background gene set files for GSEA enrichment analysis, c2.cp.kegg.v7.4.symbols.gmt c5.go.v7.4.symbols.gmt which were downloaded from http://www.gsea-msigdb.org/gsea/downloads_archive.jsp Verification box diagram GSE54236 data gSE54236_series_matrix.txt. gz was downloaded from the GEO official website (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54236/matrix/), the platform file is GPL6480-9577. TXT

## Contact us

If you have any questions, please send us an email at sl123456@qlu.edu.cn

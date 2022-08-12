# The pipeline of Hepatocellular carcinoma bioinformatics analysis
In the study “An Adenylate Uridylate (AU)-Rich Element Genes Related Prognosis Model for Hepatocellular Carcinoma based on Bioinformatics”, the RNA-seq data and clinical information of HCC patients were acquired from The Cancer Genome Atlas (TCGA) database and Gene Expression Omnibus (GEO) database. A total of 189 differentially expressed Adenylate Uridylate (AU)-Rich Element genes (DE-AREGs) between normal and HCC samples were identified in HCC. Prognostic related DE-AREGs were selected by univariate Cox regression analysis and least absolute shrinkage and selection operator (LASSO) cox analysis, including CENPA, TXNRD1, RABIF, UGT2B15, and SERPINE1. Furthermore, multivariate Cox regression analysis was used to construct an Adenylate Uridylate (AU)-Rich Element genes (AREGs) related signature. Moreover, Kaplan-Meier curves and the receiver operating characteristic (ROC) curves revealed that the AREGs-related signature could effectively predict the prognosis of HCC patients. Additionally, T stage and risk score were screened as independent prognostic factors and C-index of 1, 3 and 5 years indicating that the nomogram performed well. Functional analysis revealed that 2242 Gene Ontology (GO) and 72 Kyoto Encyclopedia of Genes and Genomes (KEGG) pathways were enriched in the high-risk group, but the low-risk group was enriched for only 2 GO terms and 13 KEGG pathways. Immune-related analyses indicated that T cell and B cell receptor abundance, microvascular endothelial cells (MVE), lymphatic endothelial cells (lye), pericytes, and stromal cells and the six immune checkpoints were significantly different between the high- and low-risk groups. In conclusion, an AREGs-related signature based on five DE-AREGs was constructed and could act as a prognostic indicator of HCC patients. Here, we will provide the steps of these analyses and some of the key document information used in these analyses.
## Data downloading and processing
### TCGA data downloading
The TCGA – LIHC dataset was retrieved tfrom the GDC official website (https://portal.gdc.cancer.gov/repository). Firstly, the FPKM file of gene expression profile was added to cart, then the cart files and metadata files were downloaded directly. The metadata file downloaded at that time can be found in attachment named as Metadata.cart.2021-12-15.json, the downloaded data is then consolidated into an expression matrix using putFilestoonedir.pl and mrna_merge.pl scripts。
``` perl
#Switch to the directory where the expression matrix is located
perl path/to/putFilestoonedir.pl #Put all gz files in one directory 
# Go to the directory obtained in the previous step
perl path/to/mrna_merge.pl path/to/Metadata.cart.2021-12-15.json #Merge the expression data of a single sample into a matrix
After that, through ensemblToSymbol.pl script in combination with Homo_sapiens.GRCh38.104.chr.gtf (http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/) annotation files, the gene ID in the matrix was converted to mRNA.symbol.txt finally.

``` perl 
path/to/ensemblToSymbol.pl path/to/Homo_sapiens.GRCh38.104.chr.gtf path/to/ensemblmatrix.txt path/to/symbolmatrix.txt #convert gene ID
Similarly, the clinical data of TCGA-LIHC dataset in GDC official website were added to CART, and the cart file was downloaded directly, and then the clinical information was gotten via getClinical.pl.
``` perl 
#Switch to the directory where the clinical information is located
perl path/to/getClinical.pl  #Extract clinical information
### GEO data downloading
GSE14520-gPL3921_series_matrix.txt.gz data obtained using GPL3921 platform were downloaded from the GEO official website (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/), the expression matrix was extracted and the matrix file GSE14520_probeExpr.txt was obtained. At the same time, the GPL3921 platform annotation file gPL3921-25447.txt was downloaded, and the probe annotation file gPL3921-25447 probe-symble.TXT was obtained. Moreover, the supplementary file gse14520_extra_supplements.txt.gz of the GSE14520 dataset was downloaded and the clinical information gse14520_clinical.txt was obtained.
## The above files are the initial input files of the data analysis script code.R. All subsequent analyses were performed using the above data.
Install R and RSrudio software, then open RSrudio software and execute the code in code.R by line.
## Other key information
Background gene set files (c2.cp.kegg.v7.4.symbols.gmt and c5.go.v7.4.symbols.gmt ) for Gene Set enrichment Analysis  were downloaded from MSigDB database (http://www.gsea-msigdb.org/gsea/downloads_archive.jsp). Moreover, gSE54236_series_matrix.txt. gz of the verification dataset (GSE54236 dataset ) was downloaded from the GEO official website (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54236/matrix/),and the platform file of the GSE54236 dataset is GPL6480-9577.txt
## Contact us
If you have any questions, please send us an email at sl123456@qlu.edu.cn

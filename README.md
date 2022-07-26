# the-data-of-Hepatocellular-Carcinoma-Bioinformatics
TCGA data sources, retrieve the TCGA – LIHC from GDC official website https://portal.gdc.cancer.gov/repository, and add the FPKM file of gene expression profile to cart, then download the cart files and metadata files directly, the metadata file downloaded at that time can be found in attachment Metadata.cart.2021-12-15.json, the downloaded data is then consolidated into an expression matrix using putFilestoonedir.pl and mrna_merge.pl scripts, after that, through ensemblToSymbol. Pl script in combination with Homo_sapiens.GRCh38.104.chr.gtf (http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/) annotation files, the gene ID in the matrix were converted to obtain mRNA.symbol.txt finally. Similarly, search TCGA-LIHC from GDC official website, add clinical data to CART, download the cart file directly, and then get the clinical information via getClinical.pl.

GEO data sources, download GSE14520-gPL3921_series_matrix.txt.gz data that use GPL3921 platform from GEO official website (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/), extract the expression matrix and obtain the matrix file GSE14520_probeExpr. TXT, at the same time, download the GPL3921 platform annotation file gPL3921-25447.txt, then the probe annotation file gPL3921-25447 probe-symble. TXT was obtained. At the same time, download the supplementary file gse14520_extra_supplements.txt. gz of this data set and extract clinical information to obtain gse14520_clinical.txt.

The above is the initial input file of the data analysis script code. R, and all subsequent analysis is based on the above data.

Background gene set files for GSEA enrichment analysis,
c2.cp.kegg.v7.4.symbols.gmt
c5.go.v7.4.symbols.gmt
which were downloaded from http://www.gsea-msigdb.org/gsea/downloads_archive.jsp
Verification box diagram GSE54236 data gSE54236_series_matrix.txt. gz was downloaded from the GEO official website (https://ftp.ncbi.nlm.nih.gov/geo/series/GSE54nnn/GSE54236/matrix/), the platform file is GPL6480-9577. TXT.

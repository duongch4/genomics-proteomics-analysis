# genomics-proteomics-analysis

## Report
My summary report can be found at "report.pdf", showing results and discussion.

## No data provided
I can only provide the R scripts that I used for the analyses, without the actual data. So they cannot be run, and for reading purpose only.
The detailed comments in each script show how to set up directories and libraries.

## My typical work flow for each of the genomic and proteomic analysis
1. Set up work environment and libraries/dependencies
2. Extract features from input data, inspect them using graphs, detect any data anomalies (if any) for quality control purpose, and prepare some solutions for those
3. Train an appropriate model from the features (files with "model" keyword), and be constantly aware of the anomalies (if any)
4. Do analysis on the differential gene expression (dge) result from step 3 (files with "furtherAnalysis" keyword):
    1. Provide visualisations on the result: word cloud, heatmap, volcano plot
    2. Network level: protein-protein interaction map, and word clouds for: GO terms, KEGG terms, protein family (Pfam), and InterPro
5. Do Gene Set Enrichment Analysis (files with "GSEA" keyword) on the dge result from step 3, and generate raw GSEA data
6. Use the raw GSEA data to provide the gene set perspective on the network level (files with "GKPI" keyword); provide network maps and word clouds for GO terms, KEGG terms, Pfam, and InterPro
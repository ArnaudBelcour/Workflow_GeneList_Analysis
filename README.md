# Workflow for Genes List Analysis

Workflow for gene enrichment analysis.

Actually, the first scripts are for GO terms enrichment analysis.

primaryFileManagement.py translates GO label from Blast2Go files into GO number.

hypergeometricTestOnGOTerms.py computes an hypergeometric test for GO terms, it also calculates different multiple tests corrections (Bonferroni, Holm, Sidak and Benjamini & Hochberg).

goTermExtractionUniprot.py queries Uniprot to obtain GO terms associated with Uniprot ID (actually the SPARQL queries don't work, only the http requests work).

This workflow works with three directories (inputFiles, temporaryFiles and outputFiles) :
	-inputFiles must have three files : queryResults.csv (a csv file resulting from queriyng on the Gene Ontology to have all the GO terms with their GO labels (it will be automated in the future), GOTermsPlasmoGenome.tsv (which contains all the GO terms from the genome of your species, it will also be automated) and your data (in .csv)).
	-temporaryFiles : will contains file used modified during the script.
	-outputFiles : the results of the analysis in tsv.
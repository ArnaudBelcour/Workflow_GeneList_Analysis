# Workflow for Genes List Analysis

Workflow for gene enrichment analysis (at this time, it only does a GO terms Enrichment analysis).

## Installation

The scripts need :

* [Python 3.0](https://www.python.org/downloads/)

Some python modules :

* [pandas](http://pandas.pydata.org/)

* [scipy](https://pypi.python.org/pypi/scipy)


# GO terms Enrichment Analysis

Actually, the first scripts are for GO terms enrichment analysis.

primaryFileManagement.py translates GO label from Blast2Go files into GO number.

goEnrichmentAnalysis.py computes an hypergeometric test for GO terms, it also calculates different multiple tests corrections (Bonferroni, Holm, Sidak and Benjamini & Hochberg).

goTermExtractionUniprot.py queries Uniprot to obtain GO terms associated with Uniprot ID (actually the SPARQL queries don't work, only the http requests work).

This workflow works with three directories (inputFiles, temporaryFiles and outputFiles) :

* inputFiles must have three files : queryResults.csv (a csv file resulting from queriyng on the Gene Ontology to have all the GO terms with their GO labels (it will be automated in the future), GOTermsPlasmoGenome.tsv (which contains all the GO terms from the genome of your species, it will also be automated) and your data (in .csv).

* temporaryFiles : will contain files used during the script, it will be created during the analysis.

* outputFiles : the results of the analysis in tsv, it will be created during the analysis.

Test used :

* Hypergeometric test to compare the distribution of GO terms in your list and in the complete organism.

And multiple testing corrections :

* [Bonferroni Correction](http://www.jstor.org/stable/2282330?seq=1#page_scan_tab_contents)

* [Sidak](https://www.jstor.org/stable/2283989?seq=1#page_scan_tab_contents)

* [Holm](http://www.jstor.org/stable/4615733?seq=1#page_scan_tab_contents)

* [Benjamini & Hochberg](https://www.jstor.org/stable/2346101?seq=1#page_scan_tab_contents)


# Launch the analysis

## Launch GO terms Enrichment Analysis

First, create a directory named "inputFiles". Put your data on it (see part "GO terms Enrichment Analysis" for more details on the needed files).
To launch the analysis, just launch the "goEnrichmentAnalysis.py" script with python.

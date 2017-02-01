# Workflow for Genes List Analysis

Workflow for gene enrichment analysis (at this time, it only does a GO terms Enrichment analysis).

## Installation

The scripts need :

* [Python](https://www.python.org/downloads/)

* [R](https://cran.r-project.org/index.html)

These python modules :

* [pandas](http://pandas.pydata.org/)

* [requests](http://docs.python-requests.org/en/master/)

* [scipy](https://pypi.python.org/pypi/scipy)

This R package :

* [KEGGREST](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html)

# Enrichment Analysis

This workflow performs enrichment analysis (especially on GO terms).

primaryFileManagement.py translates GO label from Blast2Go files into GO number.

enrichmentAnalysis.py is dividied in tow class. The first class ("EnrichmentAnalysis") is the basic method, which computes an hypergeometric test for variables (now it works for GO terms, in the future it will work for pathway) and calculates different multiple tests corrections (Bonferroni, Holm, Sidak, Benjamini & Hochberg and SGoF). The second class ("GOEnrichmentAnalysis") inherits from "EnrichmentAnalysis" and overrides a function to add GO label to the results.

goTermExtractionUniprot.py queries Uniprot to obtain GO terms associated with Uniprot ID.

This workflow works with three directories (inputFiles, temporaryFiles and outputFiles) :

* inputFiles must have three files : queryResults.csv (a csv file resulting from queriyng on the Gene Ontology to have all the GO terms with their GO labels (it will be automated in the future), GOTermsPlasmoGenome.tsv (which contains all the GO terms from the genome of your species, it will also be automated) and your data (in .csv).

* temporaryFiles : will contain files used during the script, it will be created during the analysis.

* outputFiles : the results of the analysis in tsv, it will be created during the analysis.

Test used :

* Hypergeometric test to compare the distribution of GO terms in your list and in the complete organism.

And multiple testing corrections :

* [Bonferroni Correction](http://www.jstor.org/stable/2282330?seq=1#page_scan_tab_contents) Dunn, Olive Jean. “Multiple Comparisons Among Means.” Journal of the American Statistical Association, vol. 56, no. 293, 1961, pp. 52–64.

* [Sidak](https://www.jstor.org/stable/2283989?seq=1#page_scan_tab_contents) Sidak, Zbynek. “Rectangular Confidence Regions for the Means of Multivariate Normal Distributions.” Journal of the American Statistical Association, vol. 62, no. 318, 1967, pp. 626–633.

* [Holm](http://www.jstor.org/stable/4615733?seq=1#page_scan_tab_contents) Holm, Sture. “A Simple Sequentially Rejective Multiple Test Procedure.” Scandinavian Journal of Statistics, vol. 6, no. 2, 1979, pp. 65–70.

* [Benjamini & Hochberg](https://www.jstor.org/stable/2346101?seq=1#page_scan_tab_contents) Benjamini, Yoav, and Yosef Hochberg. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” Journal of the Royal Statistical Society. Series B (Methodological), vol. 57, no. 1, 1995, pp. 289–300.

* [SGoF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2719628/) Carvajal-Rodríguez, Antonio, Jacobo de Uña-Alvarez, and Emilio Rolán-Alvarez. “A New Multitest Correction (SGoF) That Increases Its Statistical Power When Increasing the Number of Tests.” BMC Bioinformatics 10 (2009): 209.

# Network Enrichment Analysis

enzymeToPathway.R translates Enzyme Commission number (ec) into pathway using KEGGREST, it creates a tsv file (in temporaryFiles directory) containing ec number associated with the pathway name and the pathway ID on KEGG.

# Launch the analysis

## Launch GO terms Enrichment Analysis

First, create a directory named "inputFiles". Put your data on it (see part "GO terms Enrichment Analysis" for more details on the needed files).
To launch the analysis, just launch the "enrichmentAnalysis.py" script with python.

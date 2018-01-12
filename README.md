# Workflow for Genes List Analysis

**WARNING:**

**This workflow was a test during an internship. Currently it is not working and I will not update it, so don't use it.**



[![Build Status](https://travis-ci.org/ArnaudBelcour/Workflow_GeneList_Analysis.svg?branch=master)](https://travis-ci.org/ArnaudBelcour/Workflow_GeneList_Analysis)
[![Coverage Status](https://coveralls.io/repos/github/ArnaudBelcour/Workflow_GeneList_Analysis/badge.svg)](https://coveralls.io/github/ArnaudBelcour/Workflow_GeneList_Analysis)

Workflow to add information to a genome and to perfom a Singular Enrichment Analysis.

## Installation

The scripts need :

* [Python](https://www.python.org/downloads/)

* [R](https://cran.r-project.org/index.html)

These python modules :

* [numpy](https://pypi.python.org/pypi/numpy)

* [pandas >=0.19.2](http://pandas.pydata.org/)

* [requests](http://docs.python-requests.org/en/master/)

* [scipy](https://pypi.python.org/pypi/scipy)

* [SPARQLWrapper](https://rdflib.github.io/sparqlwrapper/)

* [tqdm](https://pypi.python.org/pypi/tqdm)

This R package :

* [KEGGREST](https://bioconductor.org/packages/release/bioc/html/KEGGREST.html)

# Enrichment Analysis

This workflow performs a Singular Enrichment Analysis. This analysis takes a list of genes (for example differentially expressed genes) and compute an enrichment term for each annotation term in this list.
For a better definition, read the [article writed by Huang et al. (2009)](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gkn923).

The idea behind this workflow is to make a tool to analyze genome and with few datas. Also, by using class it allows to use separately some of its components (like enrichmentAnalysis class).

fileManagement.py manages the files to create counting files used by enrichmentAnalysis.py.

enrichmentAnalysis.py is dividied in tow class. The first class ("EnrichmentAnalysis") is the basic method, which computes an hypergeometric test for variables (now it works for GO terms, in the future it will work for pathway) and calculates different multiple tests corrections (Bonferroni, Holm, Sidak, Benjamini & Hochberg and SGoF). The second class ("GOEnrichmentAnalysis") inherits from "EnrichmentAnalysis" and overrides a function to add GO label to the results.

workflow_manager.py is the main script.

This workflow works with three directories (inputFiles, temporaryFiles and outputFiles) :

* inputFiles must have three files : queryResults.csv (a csv file resulting from queriyng on the Gene Ontology to have all the GO terms with their GO labels (it will be automated in the future), GOTermsPlasmoGenome.tsv (which contains all the GO terms from the genome of your species, it will also be automated) and your data (in .csv).

* temporaryFiles : will contain files used during the script, it will be created during the analysis. It contains datas extracted from the external databases.

* outputFiles : the results of the analysis in tsv, it will be created during the analysis.

Test used :

* Hypergeometric test to compare the distribution of GO terms in your list and in the complete organism.

* Normal approximation when using big numbers.

And multiple testing corrections :

* [Bonferroni Correction](http://www.jstor.org/stable/2282330?seq=1#page_scan_tab_contents) Dunn, Olive Jean. “Multiple Comparisons Among Means.” Journal of the American Statistical Association, vol. 56, no. 293, 1961, pp. 52–64.

* [Sidak](https://www.jstor.org/stable/2283989?seq=1#page_scan_tab_contents) Sidak, Zbynek. “Rectangular Confidence Regions for the Means of Multivariate Normal Distributions.” Journal of the American Statistical Association, vol. 62, no. 318, 1967, pp. 626–633.

* [Holm](http://www.jstor.org/stable/4615733?seq=1#page_scan_tab_contents) Holm, Sture. “A Simple Sequentially Rejective Multiple Test Procedure.” Scandinavian Journal of Statistics, vol. 6, no. 2, 1979, pp. 65–70.

* [Benjamini & Hochberg](https://www.jstor.org/stable/2346101?seq=1#page_scan_tab_contents) Benjamini, Yoav, and Yosef Hochberg. “Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing.” Journal of the Royal Statistical Society. Series B (Methodological), vol. 57, no. 1, 1995, pp. 289–300.

* [Benjamini & Yekutieli](http://www.jstor.org/stable/2674075?seq=1#page_scan_tab_contents) Benjamini, Yoav, and Daniel Yekutieli. “The Control of the False Discovery Rate in Multiple Testing under Dependency.” The Annals of Statistics, vol. 29, no. 4, 2001, pp. 1165–1188.

* [SGoF](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2719628/) Carvajal-Rodríguez, Antonio, Jacobo de Uña-Alvarez, and Emilio Rolán-Alvarez. “A New Multitest Correction (SGoF) That Increases Its Statistical Power When Increasing the Number of Tests.” BMC Bioinformatics 10 (2009): 209.

# Genome enrichment annotation

This part of the workflow queries some databases to retrieve external information (GO terms, EC numbers, IPR, pathways from KEGG, MetaCyc and Reactome):

chebi_from_go.py translates GO into ChEBI using GO-plus.owl file.

database_mapping_from_gos.pys retrieves IPR, EC numbers, KEGG and MetaCyc pathways with GO terms.

eupathdb_pathway_extraction.py translates EC numbers and ChEBI into KEGG and MetaCyc pathways.

ghost_koala_pathway_extraction.py extracts KEGG pathways from Ghost Koala results.

interpro_pathway_extraction.py retrieves MetaCyc, KEGG and Reactome pathways with IPR.

keggrest_pathway_extraction.R translates Enzyme Commission number (ec) into pathway using KEGGREST, it creates a tsv file (in temporaryFiles directory) containing ec number associated with the pathway name and the pathway ID on KEGG.

panther_pathway_mapping_uniprot.py retrieves pathways from Uniprot ID.

reactome_pathway_extraction.py retrieves Reactome pathway with GO terms, EC numbers, Interpro IDs and ChEBI).

sparql_query_reactome_pathway_name.py extracts pathway name and pathway ID of Reactome pathway.

uniprot_retrieval_data.py retrieves data (GO terms, IPR and EC numbers) using Ensembl ID from Blast results.

# Launch the analysis

## Launch GO terms Enrichment Analysis

First, create a directory named "inputFiles". Put your data on it (see part "GO terms Enrichment Analysis" for more details on the needed files).
To launch the analysis, just launch the "workflow_manager.py" script with python.

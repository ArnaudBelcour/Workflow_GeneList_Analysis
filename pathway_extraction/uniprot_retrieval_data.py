#!/usr/bin/env python3

import math
import pandas as pa
import six

from SPARQLWrapper import SPARQLWrapper, JSON

from . import *

def extract_information_from_uniprot(results_dataframe):
    '''
        Requests the SPARQL endpoint of Uniprot to retrieve (from Ensembl transcrit ID) GO terms, interpro, pfam/supfam and prosites.
        The file taken as input file contains each gene associated with the result of a blast (that's the thing with 'hypothetical protein').
    '''
    results_dataframe['Blast'] = results_dataframe['Blast'].str[len('CEP03957.1hypothetical protein '):]
    results_dataframe['Blast'] = results_dataframe['Blast'].str.replace(", partial", "")

    results_dataframe.set_index("Gene_Name", inplace = True)

    for gene, row in results_dataframe.iterrows():
        transcript = 'ensembl:' + row['Blast']
        gos_found = []
        sparql = SPARQLWrapper('http://beta.sparql.uniprot.org/sparql')
        sparql.setQuery("""
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX up:<http://purl.uniprot.org/core/>
        PREFIX ensembl:<http://rdf.ebi.ac.uk/resource/ensembl/>

        SELECT DISTINCT ?go
        WHERE
        {
            ?transcrit up:transcribedFrom  ?ensemblName.
            ?protein rdfs:seeAlso ?transcrit .
            ?protein up:classifiedWith ?go .
            FILTER (regex(str(?go), "GO")) .
            VALUES ?ensemblName {""" + transcript + """}
        }
        """)

        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

        for result in results["results"]["bindings"]:
            gos_found.append(result["go"]["value"][31:])

        sparql.setQuery("""
        PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#> 
        PREFIX up:<http://purl.uniprot.org/core/>
        PREFIX ensembl:<http://rdf.ebi.ac.uk/resource/ensembl/>

        SELECT DISTINCT ?data
        WHERE
        {
            ?transcrit up:transcribedFrom  ?ensemblName.
            ?protein rdfs:seeAlso ?transcrit .
            ?protein rdfs:seeAlso ?data .
            VALUES ?ensemblName {""" + transcript + """}
        }
        """)

        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()
        datas_found = []
        print(results["results"]["bindings"])
        for result in results["results"]["bindings"]:
            datas_found.append(result["data"]["value"][len('http://purl.uniprot.org/'):])

        interpros = []
        supfams = []
        pfams = []
        prosites = []

        for data in datas_found:
            if 'interpro' in data:
                data = data[len('interpro/'):]
                interpros.append(data)
            if 'supfam' in data:
                data = data[len('supfam/'):]
                supfams.append(data)
            if 'pfam' in data and 'supfam' not in data:
                data = data[len('pfam/'):]
                pfams.append(data)
            if 'prosite' in data:
                data = data[len('prosite/'):]
                prosites.append(data)

        if type(row['GOs']) is float:
            results_dataframe.set_value(gene, 'GOs', str(gos_found))
        if row['InterProScan'] == 'no IPS match':
            results_dataframe.set_value(gene, 'InterProScan', str(interpros))

        #results_dataframe.set_value(gene, 'supFams', str(supfams))
        #results_dataframe.set_value(gene, 'pfams', str(pfams))
        #results_dataframe.set_value(gene, 'prosites', str(prosites))

    return results_dataframe

#!/usr/bin/env python3

import csv
import pandas as pa
from ast import literal_eval
from SPARQLWrapper import SPARQLWrapper, JSON

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
output_directory = 'outputFiles/'

def extract_information_from_uniprot():
    '''
        Requests the SPARQL endpoint of Uniprot to retrieve (from Ensembl transcrit ID) GO terms, interpro, pfam/supfam and prosites.
        The file taken as input file contains each gene associated with the result of a blast (that's the thing with 'hypothetical protein').
    '''
    df_genome = pa.read_excel('inputFiles/Annotation_blast2go_PROT_eH.xls', sep = "\t")
    df_genome['Blast_sur_e3'] = df_genome['Blast_sur_e3'].str[len('CEP03957.1hypothetical protein '):]
    df_genome['Blast_sur_e3'] = df_genome['Blast_sur_e3'].str.replace(", partial", "")
    df_genome = df_genome[['SeqName_eH', 'Blast_sur_e3']]
    df_genome.set_index("SeqName_eH")
    df_genome = df_genome[:50]
    for gene, row in df_genome.iterrows():
        transcript = 'ensembl:' + row['Blast_sur_e3']
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

        df_genome.set_value(gene, 'Interpros', str(interpros))
        df_genome.set_value(gene, 'GOs', str(gos_found))
        df_genome.set_value(gene, 'supFams', str(supfams))
        df_genome.set_value(gene, 'pfams', str(pfams))
        df_genome.set_value(gene, 'prosites', str(prosites))

    df_genome.to_csv(temporary_directory + 'annotations_from_uniprot.tsv', sep = "\t")
extract_information_from_uniprot()

#!/usr/bin/env python3

import csv

from SPARQLWrapper import SPARQLWrapper, JSON
from tqdm import tqdm

from . import *

def go_to_chebi():
    '''
        Requests a SPARQL endpoint in which the go-plus.owl file (from the Gene Ontology site) has been load into it.
        The query looks for all the ChEBI ID associated with a GO term.
        Then the results are writed into a csv file.
    '''
    sparql = SPARQLWrapper('http://localhost:3030/go/query')
    sparql.setQuery("""
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>

    SELECT ?go ?chebi
    WHERE {
        ?go rdfs:subClassOf  ?node .
        ?node owl:someValuesFrom  ?chebi .
        FILTER (regex(str(?chebi), "CHEBI"))
    }
    """)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    csvfile = open(temporary_directory_database + "go_chebi_mapping.tsv", "w", newline="")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['GOs', 'ChEBI'])

    for result in tqdm(results["results"]["bindings"]):
        go = result["go"]["value"][len('http://purl.obolibrary.org/obo/'):].replace("_", ":")
        chebi = result["chebi"]["value"][len('http://purl.obolibrary.org/obo/'):].replace("_", ":")
        writer.writerow([go, chebi])

    csvfile.close()


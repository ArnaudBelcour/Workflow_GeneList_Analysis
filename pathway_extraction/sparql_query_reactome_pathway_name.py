#!/usr/bin/env python3

import csv

from SPARQLWrapper import SPARQLWrapper, JSON
from tqdm import *

from . import *

def sparql_query(sparql_endpoint, query, output_file):
    '''
        Requests Reactome SPARQL Endpoint to extract correspondence between Pathway ID and Pathway Name.
        Then write it into a tsv file.
    '''
    sparql = SPARQLWrapper(sparql_endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)

    results = sparql.query().convert()

    column_names = []
    for result in results["results"]["bindings"]:
        for key in result:
            if key not in column_names:
                column_names.append(key)

    csvfile = open(temporary_directory_database + output_file + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")

    writer.writerow((column_names))

    for result in tqdm(results["results"]["bindings"]):
        writer.writerow((result["pathway_REACT"]["value"][len('http://identifiers.org/reactome/'):],
                         result["pathway_R"]["value"], result["pathway_Name"]["value"]))

    csvfile.close()

def main():
    query_file = open('sparql_queries/reactome_pathway_query.sparql', 'r')
    query = query_file.read()
    query_file.close()

    sparql_query("https://www.ebi.ac.uk/rdf/services/reactome/sparql", query, "pathwayReactomeIDToPathwayName")



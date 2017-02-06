#!/usr/bin/env python3

import csv
from SPARQLWrapper import SPARQLWrapper, JSON

temporary_directory_database = '../temporaryFiles/databases/'

def sparql_query(sparql_endpoint, query, output_file):
    sparql = SPARQLWrapper(sparql_endpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)

    results = sparql.query().convert()

    column_names = []
    for result in results["results"]["bindings"]:
        for key in result:
            if key not in column_names:
                column_names.append(key.encode("utf-8"))

    csvfile = open(temporary_directory_database + output_file + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")

    writer.writerow((column_names))

    for result in results["results"]["bindings"]:
        writer.writerow((result["pathwayID"]["value"][len('http://identifiers.org/reactome/'):].encode("utf-8"), result["pathwayName"]["value"].encode("utf-8")))

    csvfile.close()

def main():
    query_file = open('../sparqlQueries/reactomePathwayQuery.sparql', 'r')
    query = query_file.read()
    query_file.close()

    sparql_query("https://www.ebi.ac.uk/rdf/services/reactome/sparql", query, "pathwayReactomeIDToPathwayName")

main()

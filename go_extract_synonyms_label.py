#!/usr/bin/env python3

import csv
from SPARQLWrapper import SPARQLWrapper, JSON

input_directory = "inputFiles/"

def go_term_synonyms(query):
    '''
        Requests a SPARQL Endpoint (which contains the go.owl file from the Gene Ontology) to extract association between GO terms, their labels and their synonyms.
        Writes it into a tsv file.
    '''
    csvfile = open(input_directory + "query_results.tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["subject", "label", "NarrowSynonym", "BroadSynonym", "RelatedSynonym"])

    sparql = SPARQLWrapper("http://localhost:3030/datanase/query")
    sparql.setQuery(query)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        subject = result["subject"]["value"][31:]
        label = result["label"]["value"]
        if "NarrowSynonym" in result:
            narrow_synonym = result["NarrowSynonym"]["value"]
        else:
            narrow_synonym = ""
        if "BroadSynonym" in result:
            broad_synonym = result["BroadSynonym"]["value"]
        else:
            broad_synonym = ""
        if "RelatedSynonym" in result:
            related_synonym = result["RelatedSynonym"]["value"]
        else:
            related_synonym = ""
        writer.writerow([subject, label, narrow_synonym, broad_synonym, related_synonym])

    csvfile.close()

query_file = open('sparqlQueries/go_label_with_synonym.sparql', 'r')
query = query_file.read()
query_file.close()

go_term_synonyms(query)

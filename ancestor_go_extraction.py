#!/usr/bin/env python3

import pandas as pa

from functools import lru_cache
from SPARQLWrapper import SPARQLWrapper, JSON

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
output_directory = 'outputFiles/'

@lru_cache(maxsize = 12400)
def go_term_ancestor(go):
    '''
        Requests a SPARQL Endpoint (which contains the go.owl file from the Gene Ontology) to retrieve GO ancestors of a specific go term.
    '''
    go_ancestors = []
    sparql = SPARQLWrapper("http://localhost:3030/datanase/query")
    sparql.setQuery("""
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX uniprot: <http://bio2rdf.org/uniprot:>
    PREFIX go: <http://purl.obolibrary.org/obo/GO_>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>

    SELECT DISTINCT ?goAnc ?goAncLabel
    WHERE {
      go:""" + go[3:] + """ (rdfs:subClassOf|(rdfs:subClassOf/owl:someValuesFrom))*
    ?goAnc .
      OPTIONAL { ?goAnc rdfs:label ?goAncLabel .}
        FILTER ((str(?goAnc) != "b")) .
    }
    """)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        go_ancestors.append(result["goAnc"]["value"][31:])

    return go_ancestors

def union_go_and_their_ancestor(gos):
    '''
        Takes a list of GO terms corresponding to the value of a dataframe column as input (here the GO terms associated with a gene).
        Uses go_term_ancestor() function to retrieve all the ancestors of the go term.
        Makes the union of the results and return the result to replace the former list in the dataframe.
    '''
    go_ancestors_for_go_lists = []

    for go in gos:
        go_ancestors = go_term_ancestor(go)
        go_ancestors_for_go_lists.append(go_ancestors)

    go_list_for_entity = list(set().union(*go_ancestors_for_go_lists))

    return go_list_for_entity

#!/usr/bin/env python3
import csv
import pandas as pa

from SPARQLWrapper import SPARQLWrapper, JSON

output_directory = 'outputFiles/'

def insert_delete_pval(prefix, go_pval_rdf, insert_or_delete):
    '''
        Insert datas into a Fuseki SPARQL Endpoint.
    '''

    sparql = SPARQLWrapper("http://localhost:3030/go/update")
    sparql.method = 'POST'
    sparql.setQuery(prefix + "\n" + insert_or_delete + """ DATA
    {
    """ + go_pval_rdf + """
    }
    """)
    sparql.query()

def go_term_ancestor(go):
    '''
        Requests a SPARQL Endpoint (which contains the go.owl file from the Gene Ontology) to retrieve GO ancestors of a specific go term.
    '''
    go_ancestors = []
    sparql = SPARQLWrapper("http://localhost:3030/go/query")
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
        FILTER contains(str(?goAnc), "GO_")
    }
    """)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        go_ancestors.append(result["goAnc"]["value"][31:].replace("_",":"))
    go_ancestors.remove(go)
    return go_ancestors

def middle_selection(under_over):
    '''
        Remove GO terms which have a more significant descendant and GO terms which have a more significant ancestor.
    '''
    df_over = pa.read_csv(output_directory + 'pValuesOfGOs_' + under_over +'.tsv', sep='\t', header=1)
    go_pval_rdf = ""
    prefix = """PREFIX od: <http://test/od/>
                PREFIX go: <http://purl.obolibrary.org/obo/>"""
    for index, row in df_over.iterrows():
        if row['pValueBenjaminiHochberg'] < 0.05:
            go_pval_rdf += "go:" + str(row['GOs']) + " " + "od:hasPValueHypergeometric" + " " + str(row['pValueBenjaminiHochberg']) + " .\n"

    go_pval_rdf = go_pval_rdf.replace('GO:', 'GO_')

    insert_delete_pval(prefix, go_pval_rdf, "INSERT")

    sparql = SPARQLWrapper("http://localhost:3030/go/query")
    sparql.setQuery("""
    PREFIX od: <http://test/od/>
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>
    PREFIX go: <http://purl.obolibrary.org/obo/>
    SELECT DISTINCT ?goID ?pValue ?goLabel
    WHERE {
        ?goID od:hasPValueHypergeometric ?pValue .
        ?goID rdfs:label ?goLabel .
        # remove all the terms that have a more significant descendant
        FILTER NOT EXISTS {
            ?goDescendant ( rdfs:subClassOf|(rdfs:subClassOf/owl:someValuesFrom))* ?goID .
            ?goDescendant od:hasPValueHypergeometric ?descendentPValue .
            FILTER (?descendentPValue < ?pValue) .
            }
        # remove all the terms that have a more significant ancestor
        FILTER NOT EXISTS {
            ?goID ( rdfs:subClassOf|(rdfs:subClassOf/owl:someValuesFrom))* ?goAncestor .
            ?goAncestor od:hasPValueHypergeometric ?ancestorPValue .
            FILTER (?ancestorPValue < ?pValue) .
            }
        }
        ORDER BY ASC(?pValue)
    """)
    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()
    csvfile = open(output_directory + "result_go_middle_" + under_over + ".tsv", "w", newline="")
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['goID', 'pValue', 'goLabel'])

    for result in results["results"]["bindings"]:
        writer.writerow([result["goID"]["value"][len('http://purl.obolibrary.org/obo/'):], 
                         result["pValue"]["value"], result["goLabel"]["value"]])
    csvfile.close()

    insert_delete_pval(prefix, go_pval_rdf, "DELETE")
 
def specific_selection(under_over):
    '''
        Remove all the ancestors of the lowest GO terms and keep these lowest GO terms.
    '''
    df_over = pa.read_csv(output_directory + 'pValuesOfGOs_' + under_over +'.tsv', sep='\t', header=1)
    significatives_gos = df_over[df_over['pValueBenjaminiHochberg'] < 0.05]['GOs'].tolist()
    go_to_delete = []

    for go in significatives_gos:
        go_with_ancestors = go_term_ancestor(go)
        go_to_delete.extend(go_with_ancestors)

    for go in go_to_delete:
        if go in significatives_gos:
            significatives_gos.remove(go)
    if significatives_gos == []:
        significatives_gos.append('')

    df_go = pa.DataFrame(significatives_gos)
    df_go.columns = [['GOs']]
    df_go.set_index('GOs', inplace=True)
    df_over.set_index('GOs', inplace=True)
    df_joined = df_go.join(df_over)
    df_joined.to_csv(output_directory + "result_go_cleaned_specific_" + under_over + ".tsv", sep="\t")

def main():
    middle_selection('over')
    specific_selection('over')
    middle_selection('under')
    specific_selection('under')

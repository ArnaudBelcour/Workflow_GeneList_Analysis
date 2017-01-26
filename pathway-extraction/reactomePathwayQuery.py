import csv
from SPARQLWrapper import SPARQLWrapper, JSON

temporaryDirectory = '../temporaryFiles/'

def sparqlQuery(sparqlEndpoint, query, outputFile):
    sparql = SPARQLWrapper(sparqlEndpoint)
    sparql.setQuery(query)
    sparql.setReturnFormat(JSON)

    results = sparql.query().convert()

    columnNames = []
    for result in results["results"]["bindings"]:
        for key in result:
            if key not in columnNames:
                columnNames.append(key.encode("utf-8"))

    csvfile = open(temporaryDirectory + outputFile + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")

    writer.writerow((columnNames))

    for result in results["results"]["bindings"]:
        writer.writerow((result["pathwayID"]["value"][len('http://identifiers.org/reactome/'):].encode("utf-8"), result["pathwayName"]["value"].encode("utf-8")))

    csvfile.close()

def main():
    queryFile = open('../sparqlQueries/reactomePathwayQuery.sparql', 'r')
    query = queryFile.read()
    queryFile.close()

    sparqlQuery("https://www.ebi.ac.uk/rdf/services/reactome/sparql", query, "pathwayReactomeIDToPathwayName")

main()

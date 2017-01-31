import csv
from SPARQLWrapper import SPARQLWrapper, JSON

def goTermSynonyms(query):
    csvfile = open("query_results.tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["subject", "label", "NarrowSynonym", "BroadSynonym", "RelatedSynonym"])

    goAncestors = []
    sparql = SPARQLWrapper("http://localhost:3030/datanase/query")
    sparql.setQuery(query)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        subject = result["subject"]["value"][31:]
        label = result["label"]["value"]
        if "NarrowSynonym" in result:
            narrowSynonym = result["NarrowSynonym"]["value"]
        else:
            narrowSynonym = ""
        if "BroadSynonym" in result:
            broadSynonym = result["BroadSynonym"]["value"]
        else:
            broadSynonym = ""
        if "RelatedSynonym" in result:
            relatedSynonym = result["RelatedSynonym"]["value"]
        else:
            relatedSynonym = ""
        writer.writerow([subject, label, narrowSynonym, broadSynonym, relatedSynonym])

    csvfile.close()

queryFile = open('sparqlQueries/goLabelWithSynonym.sparql', 'r')
query = queryFile.read()
queryFile.close()

goTermSynonyms(query)

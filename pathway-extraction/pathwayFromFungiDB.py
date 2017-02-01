import csv
import requests

temporaryDirectory = '../temporaryFiles/'

def requestFungiDB():
    metacycPathways = []
    keggPathways = []
    r = requests.get('http://fungidb.org/common/downloads/pathwayFiles//MetaCyc/')
    for row in r.text.split("\n"):
        if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
            metacycPathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

    r = requests.get('http://fungidb.org/common/downloads/pathwayFiles//KEGG/')
    for row in r.text.split("\n"):
        if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
            keggPathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

    return metacycPathways, keggPathways

def linkPathwayIDrequestFungiDB(database, pathways):
    csvfile = open(temporaryDirectory + "ecChebiToPathwayFungiDB" + database + ".tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(('pathway', 'ecChebis'))

    for pathway in pathways:
        if '.xgmml' in pathway:
            ecChebis = []
            r = requests.get('http://fungidb.org/common/downloads/pathwayFiles//' + database + '/' + pathway)
            for row in r.text.split("\n"):
                if "node label" in row:
                    ecChebis.append(row[len('  <node label="'):].split('"')[0])
            writer.writerow((pathway[:-len('.xgmml')], ecChebis))

    csvfile.close()

def main():
    metacycPathways, keggPathways = requestFungiDB()
    linkPathwayIDrequestFungiDB('KEGG', keggPathways)
    linkPathwayIDrequestFungiDB('MetaCyc', metacycPathways)

main()

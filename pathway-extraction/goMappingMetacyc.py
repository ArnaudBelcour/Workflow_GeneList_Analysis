import csv
import pandas as pa
import requests

temporaryDirectory = '../temporaryFiles/'

def httpRequestGeneOntology(url, fileName):
    r = requests.get('http://geneontology.org/external2go/metacyc2go')
    r = r.text.encode("utf-8")
    csvfile = open(temporaryDirectory + fileName + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    ar = r.split("\n")
    writer.writerow(('metacycPathway', 'goLabel', 'goNumber'))
    for r in ar:
        if r.startswith("!"):
            pass
        else:
            if 'PWY' in r:
                r = r.replace(">", ";")
                writer.writerow((r.split(" ; ")))

    csvfile.close()

    return r

def cleaningFile(fileName):
    df = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
    df = df[:-1]
    df['metacycPathway'] = df['metacycPathway'].str.replace("MetaCyc:", "")
    df['goLabel'] = df['goLabel'].str.replace("GO:", "")
    df.to_csv(temporaryDirectory + "output.tsv", sep= "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

def main():
    fileName = 'metacycGOTranslation'
    test = httpRequestGeneOntology('http://geneontology.org/external2go/metacyc2go', fileName)
    cleaningFile(fileName)

main()

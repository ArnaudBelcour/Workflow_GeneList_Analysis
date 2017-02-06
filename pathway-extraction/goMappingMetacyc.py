#!/usr/bin/env python3

import csv
import pandas as pa
import requests

temporary_directory = 'temporaryFiles/'
temporary_directory_database = '../temporaryFiles/databases/'

def http_request_gene_ontology(url, file_name):
    r = requests.get('http://geneontology.org/external2go/metacyc2go')
    r = r.text.encode("utf-8")
    csvfile = open(temporary_directory + file_name + ".tsv", "w")
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

def cleaning_file(file_name):
    df = pa.read_csv(temporary_directory_database + file_name + ".tsv", sep = "\t")
    df = df[:-1]
    df['metacycPathway'] = df['metacycPathway'].str.replace("MetaCyc:", "")
    df['goLabel'] = df['goLabel'].str.replace("GO:", "")
    df.to_csv(temporary_directory_database + file_name + ".tsv", sep= "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

def main():
    file_name = 'metacycGOTranslation'
    test = http_request_gene_ontology('http://geneontology.org/external2go/metacyc2go', file_name)
    cleaning_file(file_name)

main()

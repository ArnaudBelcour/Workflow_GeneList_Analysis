#!/usr/bin/env python3

import csv
import pandas as pa
import urllib2

temporary_directory = 'temporaryFiles/'
temporary_directory_database = '../temporaryFiles/databases/'

def http_request_gene_ontology(url, file_name):
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    page = response.read()
    page = page.split("\n")

    csvfile = open(temporary_directory + file_name + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['pathwayAccession', 'pathwayName', 'uniProtID', 'pantherSubfamilyID', 'pantherSubfamilyName'])

    for i in page[:-1]:
        column_separation = i.split("\t")
        writer.writerow([column_separation[0], column_separation[1], column_separation[4], column_separation[9], column_separation[10].replace('\r', '')])

    csvfile.close()

    return page

def cleaning_file(file_name):
    df = pa.read_csv(temporary_directory + file_name + ".tsv", sep = "\t")
    df['metacycPathway'] = df['metacycPathway'].str.replace("MetaCyc:", "")
    df['goLabel'] = df['goLabel'].str.replace("GO:", "")
    df.to_csv(temporary_directory_database + file_name + ".tsv", sep= "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

def main():
    file_name = 'uniprotPantherPathwayTranslation'
    test = http_request_gene_ontology('ftp://ftp.pantherdb.org/pathway/current_release/SequenceAssociationPathway3.4.1.txt', file_name)
    cleaning_file(file_name)

main()

#!/usr/bin/env python3

import csv
import pandas as pa

from tqdm import *
from urllib.request import urlopen

from . import *

def http_request_gene_ontology(url, file_name):
    '''
        Requests Panther ftp to obtain the mapping file between Uniprot ID/Ensembl ID and pathway.
    '''
    response = urlopen(url)
    page = response.read().decode("utf-8")
    page = page.split("\n")

    csvfile = open(temporary_directory_database + file_name + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['pathway_accession', 'pathway_name', 'uniprot_id', 'panther_subfamily_id', 'panther_subfamily_name'])

    for i in tqdm(page[:-1]):
        column_separation = i.split("\t")
        writer.writerow([column_separation[0], column_separation[1], column_separation[4], column_separation[9], column_separation[10].replace('\r', '')])

    csvfile.close()

def cleaning_file(file_name):
    df = pa.read_csv(temporary_directory_database + file_name + ".tsv", sep = "\t")
    df['uniprot_id'] = [pathway.split("=")[-1]for pathway in df['uniprot_id']]
    df.to_csv(temporary_directory_database + file_name + ".tsv", sep = "\t")

def main():
    file_name = 'uniprot_panther_pathway'
    http_request_gene_ontology('ftp://ftp.pantherdb.org/pathway/current_release/SequenceAssociationPathway3.4.1.txt', file_name)
    cleaning_file(file_name)


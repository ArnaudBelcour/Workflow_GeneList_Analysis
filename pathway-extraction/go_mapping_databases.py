#!/usr/bin/env python3

import csv
import pandas as pa
import requests

temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'

def http_request_gene_ontology(url, file_name):
    """
        Requests the Gene Ontology server to obtain mapping file between GO and Interpro, KEGG, Enzyme Code, MetaCyc.
        Rewrites each file into a correct tsv file.
    """
    response = requests.get(url)
    results = response.text
    csvfile = open(temporary_directory_database + file_name + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    results_splitted = results.split("\n")
    if file_name == 'metacyc_go_mapping':
        id_name = 'metacyc_pathway'
        id_prefix = 'MetaCyc:'
    elif file_name == 'reactome_go_mapping':
        id_name = 'reactome_pathway'
        id_prefix = 'Reactome:'
    elif file_name == 'kegg_go_mapping':
        id_name = 'kegg_pathway'
        id_prefix = 'KEGG:'
    elif file_name == 'interpro_go_mapping':
        id_name = 'interpro'
        id_prefix = 'InterPro:'
    elif file_name == 'eccode_go_mapping':
        id_name = 'ec_code'
        id_prefix = 'EC:'
    writer.writerow((id_name, 'go_label', 'GOs'))
    for row in results_splitted:
        if row.startswith("!"):
            pass
        else:
            if file_name == 'metacyc_go_mapping':
                if 'PWY' in row:
                    row = row.replace(">", ";")
                    writer.writerow((row.split(" ; ")))
            else:
                row = row.replace(">", ";")
                writer.writerow((row.split(" ; ")))

    csvfile.close()

    return id_name, id_prefix

def cleaning_file(file_name, id_name, id_prefix):
    df = pa.read_csv(temporary_directory_database + file_name + ".tsv", sep = "\t")
    df = df[:-1]
    df[id_name] = df[id_name].str.replace(id_prefix, "")
    if id_name == "interpro":
        df[id_name] = [[interpro[:9]
                                for interpro in interpros]
                            for interpros in df[id_name]]
    df['GOs'] = df['GOs'].str.replace("GO:", "GO_")
    df.to_csv(temporary_directory_database + file_name + ".tsv", sep= "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

def main():
    databases_gos_mapping = {'metacyc_go_mapping': 'http://geneontology.org/external2go/metacyc2go',
                   'reactome_go_mapping': 'http://geneontology.org/external2go/reactome2go',
                   'kegg_go_mapping': 'http://geneontology.org/external2go/kegg2go',
                   'interpro_go_mapping': 'http://geneontology.org/external2go/interpro2go',
                   'eccode_go_mapping': 'http://geneontology.org/external2go/ec2go',
                   }
    for database in databases_gos_mapping:
        id_name, id_prefix = http_request_gene_ontology(databases_gos_mapping[database], database)
        cleaning_file(database, id_name, id_prefix)

main()

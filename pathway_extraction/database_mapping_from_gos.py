#!/usr/bin/env python3

import csv
import pandas as pa

from tqdm import *

from . import *

def request_gene_ontology(url, file_name):
    """
        Requests the Gene Ontology server to obtain mapping file between GO and Interpro, KEGG, Enzyme Code, MetaCyc.
        Rewrites each file into a correct tsv file.
    """
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

    if file_name in ['metacyc_go_mapping', 'reactome_go_mapping', 'kegg_go_mapping', 'eccode_go_mapping']:
        df = pa.read_csv(url, sep=' > | ; ', skiprows=2, header=None, engine='python')
    elif file_name  == 'interpro_go_mapping':
        df = pa.read_csv(url, sep=' > | ; ', skiprows=6, header=None, engine='python')

    df.columns = [[id_name, 'go_label', 'GOs']]
    df[id_name] = df[id_name].str.replace(id_prefix, "")
    df[id_name] = df[id_name].str.strip(to_strip='+-')

    if id_name == "interpro":
        df[id_name] = [interpro[:9]
                        for interpro in df[id_name]]
    if id_name == "ec_code":
        df[id_name] = ['ec:'+ec
                        for ec in df[id_name]]

    df['go_label'] = df['go_label'].str.replace("GO:", "")
    df['go_label'] = df['go_label'].str.strip(to_strip='+-')

    df.to_csv(temporary_directory_database + file_name + ".tsv", sep='\t', index=False, header=True, quoting=csv.QUOTE_NONE)

def main():
    databases_gos_mapping = {'metacyc_go_mapping': 'http://geneontology.org/external2go/metacyc2go',
                   'reactome_go_mapping': 'http://geneontology.org/external2go/reactome2go',
                   'kegg_go_mapping': 'http://geneontology.org/external2go/kegg2go',
                   'interpro_go_mapping': 'http://geneontology.org/external2go/interpro2go',
                   'eccode_go_mapping': 'http://geneontology.org/external2go/ec2go',
                   }

    for database in tqdm(databases_gos_mapping):
        request_gene_ontology(databases_gos_mapping[database], database)

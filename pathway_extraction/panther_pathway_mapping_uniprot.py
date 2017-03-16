#!/usr/bin/env python3

import pandas as pa

from tqdm import *

from . import *

def request_panther_file():
    '''
        Requests Panther ftp to obtain the mapping file between Uniprot ID/Ensembl ID and pathway.
    '''
    url = 'ftp://ftp.pantherdb.org/pathway/current_release/SequenceAssociationPathway3.4.1.txt'
    file_name = 'uniprot_panther_pathway'

    print("\tRequesting Panther.")
    df = pa.read_csv('ftp://ftp.pantherdb.org/pathway/current_release/SequenceAssociationPathway3.4.1.txt', sep='\t', header=None)

    df.columns = [['pathway_accession', 'pathway_name', 'pathway_component_accession', 'pathway_component_name', 'uniprot_ID',
                   'protein_definition', 'confidence_code', 'evidence', 'evidence_type', 'panter_subfamily_ID',
                   'panther_subfamily_name']]
    df = df[['pathway_accession', 'pathway_name', 'uniprot_ID', 'panter_subfamily_ID', 'panther_subfamily_name']]

    df['uniprot_ID'] = [pathway.split("=")[-1]for pathway in df['uniprot_ID']]

    df.to_csv(temporary_directory_database + file_name + ".tsv", sep='\t', index=False)

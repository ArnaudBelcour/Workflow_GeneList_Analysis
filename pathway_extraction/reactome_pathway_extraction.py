#!/usr/bin/env python3

import csv
import numpy as np
import pandas as pa
import requests

from tqdm import *

from . import *

def http_request_reactome(data_id, data_name, writer, session=requests):
    '''
        Requests Reactome to retrieve pathway associated with an ID (GO terms, Reactome ID, CHEBI ID and Enzyme Code ID).
    '''

    try:
        r = session.get('http://www.reactome.org/ContentService/search/query?query=' + data_id +'&cluster=true')
        results = r.json()

        if 'code' in results:
            if results['code'] == 404:
                return

        if data_name in ["GO", "EC"]:
            data_id = data_name + ":" + data_id

        for index in range(len(results['results'])):
            for index2 in range(len(results['results'][index]['entries'])):
                reactome_id = results['results'][index]['entries'][index2]['id']
                reactome_id_specie = results['results'][index]['entries'][index2]['species']
                writer.writerow([data_id, reactome_id, reactome_id_specie])

    except Exception as e:
        error = "Errors : " + repr(e)

def file_creation(data_name, column_name, df_genome, session=requests):
    csvfile = open(temporary_directory_database + 'pathway_reactome_' + data_name + '.tsv', "w", newline="")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow([data_name, 'Id', 'specie', 'data_type'])

    datas_requests = []

    [datas_requests.extend(datas.split(",")) for datas in df_genome[column_name]]

    datas_requests = list(set().union(datas_requests))

    if data_name == "Interpro":
        print("\tPathway from Interpro extraction")
        for data in tqdm(datas_requests):
            http_request_reactome(data, data_name, writer, session)
    elif data_name == "CHEBI":
        print("\tPathway from ChEBI extraction")
        for data in tqdm(datas_requests):
            data = data.strip().replace("_", ":")
            http_request_reactome(data, data_name, writer, session)
    elif data_name in ["GO", "EC"]:
        print("\tPathway from " + data_name + " extraction")
        for data in tqdm(datas_requests):
            data = data.strip()[len(data_name + ':'):]
            http_request_reactome(data, data_name, writer, session)

    csvfile.close()

def main(file_name, session=requests):
    df_genome = pa.read_csv(temporary_directory + file_name, sep="\t")
    df_genome.replace(np.nan, '', regex=True, inplace=True)
    data_names = {"EC": "EnzymeCodes", "GO": "GOs", "Interpro": "InterProScan", "CHEBI": "ChEBI"}

    for data_name in data_names:
        file_creation(data_name, data_names[data_name], df_genome, session)

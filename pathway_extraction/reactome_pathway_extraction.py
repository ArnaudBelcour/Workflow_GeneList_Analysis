#!/usr/bin/env python3

import csv
import pandas as pa
import requests

temporary_directory_database = '../temporaryFiles/databases/'

def http_request_reactome(data_id, data_name, writer):
    '''
        Requests Reactome to retrieve pathway associated with an ID (GO terms, Reactome ID, CHEBI ID and Enzyme Code ID).
    '''

    try:
        r = requests.get('http://www.reactome.org/ContentService/search/query?query=' + data_id +'&cluster=true')
        results = r.json()
        print(results)

        if 'code' in results:
            if results['code'] == 404:
                return print('No results found')

        #r = requests.get('http://www.reactome.org/ContentService/data/pathways/low/entity/' + result_id)
        if data_name in ["REACT", "Interpro", "GO", "EC"]:
            if data_name in ["GO", "EC"]:
                data_id = data_name + "_" + data_id
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie = results['results'][0]['entries'][index]['species']

        elif data_name == "CHEBI":
            for index in range(len(results['results'])):
                for index2 in range(len(results['results'][index]['entries'])):
                    reactome_id = results['results'][index]['entries'][index2]['id']
                    reactome_id_specie = results['results'][index]['entries'][index2]['species']

        writer.writerow([data_id, reactome_id, reactome_id_specie])

    except Exception as e:
        print ("Errors : " + repr(e))

def file_creation(data_name, column_name, df_genome):
    csvfile = open(temporary_directory_database + 'pathway_reactome_' + data_name + '.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter = "\t")
    writer.writerow([data_name, 'Id', 'specie', 'data_type'])

    datas_requests = []
    if data_name in ["EC", "Interpro"]:
        [datas_requests.extend(datas.split("; ")) for datas in df_genome[column_name].dropna().tolist()]
    elif data_name in ["GO", "CHEBI", "REACT"]:
        [datas_requests.extend(datas.split(",")) for datas in df_genome[column_name]]

    datas_requests = list(set().union(datas_requests))

    if data_name == "Interpro":
        for data in datas_requests:
            data = data.strip()[:9]
            http_request_reactome(data, data_name, writer)
    elif data_name == "CHEBI":
        for data in datas_requests:
            data = data.strip().replace("_", ":")
            http_request_reactome(data, data_name, writer)
    elif data_name == "REACT":
        for data in datas_requests:
            data = data.strip()
            http_request_reactome(data, data_name, writer)
    elif data_name in ["GO", "EC"]:
        for data in datas_requests:
            data = data.strip()[len(data_name + ':'):]
            http_request_reactome(data, data_name, writer)

    csvfile.close()

def main():
    df_genome = pa.read_csv('../inputFiles/genome_with_pathway.tsv', sep = "\t")
    data_names = {"EC": "EnzymeCodes", "GO": "GOs", "CHEBI":"ChEBI", "REACT": "reactome_pathway", "Interpro": "InterProScan"}

    for data_name in data_names:
        file_creation(data_name, data_names[data_name], df_genome)
main()

#!/usr/bin/env python3

import csv
import pandas as pa
import requests
from ast import literal_eval

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
        if data_name == "GO":
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie = results['results'][0]['entries'][index]['species']
                writer.writerow([data_name + "_" + data_id, reactome_id, reactome_id_specie])
        elif data_name == "EC":
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie =results['results'][0]['entries'][index]['species']
                writer.writerow([data_name + ":" + data_id, reactome_id, reactome_id_specie])
        elif data_name == "REACT":
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie = results['results'][0]['entries'][index]['species']
                writer.writerow([data_name + "_" + data_id, reactome_id, reactome_id_specie])
        elif data_name == "CHEBI":
            for index in range(len(results['results'])):
                for index2 in range(len(results['results'][index]['entries'])):
                    reactome_id = results['results'][index]['entries'][index2]['id']
                    reactome_id_specie = results['results'][index]['entries'][index2]['species']
                    writer.writerow([data_name + "_" + data_id, reactome_id, reactome_id_specie])
        elif data_name == "Interpro":
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie = results['results'][0]['entries'][index]['species']
                writer.writerow([data_id, reactome_id, reactome_id_specie])

    except Exception as e:
        print ("Errors : " + repr(e))

def main():
    df_genome = pa.read_csv('../inputFiles/genome_with_pathway.tsv', sep = "\t")
    data_names = ["EC", "GO", "CHEBI", "REACT", "Interpro"]

    csvfile = open('../temporaryFiles/databases/enzyme_pathway_reactome_EC.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["EC", 'Id', 'specie'])
    ecs_requests = []
    [ecs_requests.extend(ecs.split("; ")) for ecs in df_genome['EnzymeCodes'].dropna().tolist()]
    ecs_requests = list(set().union(ecs_requests))
    for ec in ecs_requests:
        ec = ec.strip()[len('ec:'):]
        http_request_reactome(ec, "EC", writer)
    csvfile.close()

    csvfile = open('../temporaryFiles/databases/enzyme_pathway_reactome_GO.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["GO", 'Id', 'specie'])
    gos_requests = []
    [gos_requests.extend(literal_eval(gos)) for gos in df_genome['GOs']]
    gos_requests = list(set().union(gos_requests))
    for go in gos_requests:
        go = go.strip()[len('GO_'):]
        http_request_reactome(go, "GO", writer)
    csvfile.close()

    csvfile = open('../temporaryFiles/databases/enzyme_pathway_reactome_CHEBI.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["CHEBI", 'Id', 'specie'])
    chebis_requests = []
    [chebis_requests.extend(literal_eval(chebis)) for chebis in df_genome['ChEBI']]
    chebis_requests = list(set().union(chebis_requests))
    for chebi in chebis_requests:
        chebi = chebi.strip()[len('CHEBI:'):]
        http_request_reactome(chebi, "CHEBI", writer)
    csvfile.close()

    csvfile = open('../temporaryFiles/databases/enzyme_pathway_reactome_REACT.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["REACT", 'Id', 'specie'])
    reactomes_requests = []
    [reactomes_requests.extend(literal_eval(reactomes)) for reactomes in df_genome['reactome_pathway']]
    reactomes_requests = list(set().union(reactomes_requests))
    for reactome in reactomes_requests:
        http_request_reactome(reactome, "REACT", writer)
    csvfile.close()

    csvfile = open('../temporaryFiles/databases/enzyme_pathway_reactome_Interpro.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(["Interpro", 'Id', 'specie'])
    interpros_requests = []
    [interpros_requests.extend(interpros.split(";")) for interpros in df_genome['InterProScan']]
    interpros_requests = list(set().union(interpros_requests))
    for interpro in interpros_requests:
        interpro = interpro.strip()[:9]
        http_request_reactome(interpro, "Interpro", writer)
    csvfile.close()
main()

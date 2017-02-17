#!/usr/bin/env python3

import csv
import requests

def http_request_reactome(code, writer):
    '''
        Requests Reactome to retrieve pathway associated with an ID (GO terms, Reactome ID, CHEBI ID and Enzyme Code ID).
    '''
    try:
        r = requests.get('http://www.reactome.org/ContentService/search/query?query=' + code +'&cluster=true')
        results = r.json()
        print(results)

        if 'code' in results:
            if results['code'] == 404:
                return print('No results found')

        #r = requests.get('http://www.reactome.org/ContentService/data/pathways/low/entity/' + result_id)
        if "GO" in code:
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie = results['results'][0]['entries'][index]['species']
                writer.writerow([code, reactome_id, reactome_id_specie])
        if "." in code:
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie =results['results'][0]['entries'][index]['species']
                writer.writerow([code, reactome_id, reactome_id_specie])
        if "REACT" in code:
            for index in range(len(results['results'][0]['entries'])):
                reactome_id = results['results'][0]['entries'][index]['id']
                reactome_id_specie = results['results'][0]['entries'][index]['species']
                writer.writerow([code, reactome_id, reactome_id_specie])
        if "CHEBI" in code:
            for index in range(len(results['results'])):
                for index2 in range(len(results['results'][index]['entries'])):
                    reactome_id = results['results'][index]['entries'][index2]['id']
                    reactome_id_specie = results['results'][index]['entries'][index2]['species']
                    writer.writerow([code, reactome_id, reactome_id_specie])

    except Exception as e:
        print ("Errors : " + repr(e))

def main():
    csvfile = open('temporaryFiles/databases/enzyme_pathway_reactome.tsv', "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['code', 'Id', 'specie'])
    http_request_reactome('CHEBI:18139', writer)
    #for ec in ecs_requests:
        #ec = ec.strip()[len('ec:'):]
        #test = http_request_reactome(ec, writer)
    csvfile.close()

main()

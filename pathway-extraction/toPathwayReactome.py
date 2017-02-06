#!/usr/bin/env python3

import requests

def http_request_reactome(code):
    try:
        r = requests.get('http://www.reactome.org/ContentService/search/query?query=' + code +'&cluster=true')
        results = r.json()
        result_id = results['results'][0]['entries'][0]['stId']

        r = requests.get('http://www.reactome.org/ContentService/data/pathways/low/entity/' + result_id)
        print(r.json())
    except:
        print(r.json()['messages'][0])

    return r.text

def main():
    test = http_request_reactome('0007264')

main()

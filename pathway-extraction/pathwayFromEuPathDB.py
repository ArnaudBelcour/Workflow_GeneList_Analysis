#!/usr/bin/env python3

import csv
import requests

temporary_directory_database = '../temporaryFiles/databases/'

def request_eupathdb(db_database):
    db_database_pathways = {}
    metacyc_pathways = []
    kegg_pathways = []

    r = requests.get('http://' + db_database + '.org/common/downloads/pathwayFiles//MetaCyc/')
    for row in r.text.split("\n"):
        if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
            metacyc_pathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

    r = requests.get('http://' + db_database + '.org/common/downloads/pathwayFiles//KEGG/')
    for row in r.text.split("\n"):
        if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
            kegg_pathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

     if db_database == 'tritrypdb':
        leishcyc_pathways = []
        trypanocyc_pathways = []
        r = requests.get('http://' + db_database + '.org/common/downloads/pathwayFiles//LeishCyc/')
        for row in r.text.split("\n"):
            if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
                leishcyc_pathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

        r = requests.get('http://' + db_database + '.org/common/downloads/pathwayFiles//TrypanoCyc/')
        for row in r.text.split("\n"):
            if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
                trypanocyc_pathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

        db_database_pathways['LeishCyc'] = leishcyc_pathways
        db_database_pathways['TrypanoCyc'] = trypanocyc_pathways

    db_database_pathways['MetaCyc'] = metacyc_pathways
    db_database_pathways['KEGG'] = kegg_pathways

    return db_database_pathways

def link_pathway_id_request_eupathdb(db_database, database, pathways):
    csvfile = open(temporary_directory_database + "ecChebiToPathway_" + db_database + database + ".tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(('pathway', 'ecChebis'))

    for pathway in pathways:
        if '.xgmml' in pathway:
            ec_and_chebis = []
            r = requests.get('http://' + db_database + '.org/common/downloads/pathwayFiles//' + database + '/' + pathway)
            for row in r.text.split("\n"):
                if "node label" in row:
                    ec_and_chebis.append(row[len('  <node label="'):].split('"')[0])
            writer.writerow((pathway[:-len('.xgmml')], ec_and_chebis))

    csvfile.close()

def main():
    db_databases = ['amoebadb', 'cryptodb', 'fungidb', 'giardiadb', 'microsporidiadb', 'piroplasmadb', 'plasmodb', 'toxodb', 'trichdb', 'tritrypdb']
    for db_database in db_databases:
        db_database_pathways = request_eupathdb(db_database)
        link_pathway_id_request_eupathdb(db_database, 'KEGG', db_database_pathways['KEGG'])
        link_pathway_id_request_eupathdb(db_database, 'MetaCyc', db_database_pathways['MetaCyc'])
        if db_database == 'tritrypdb':
            link_pathway_id_request_eupathdb(db_database, 'LeishCyc', db_database_pathways['LeishCyc'])
            link_pathway_id_request_eupathdb(db_database, 'TrypanoCyc', db_database_pathways['TrypanoCyc'])

main()

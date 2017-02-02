#!/usr/bin/env python3

import csv
import requests

temporaryDirectoryDatabase = '../temporaryFiles/databases/'

def requestEuPathDB(dbDatabase):
    dbDatabasePathways = {}
    metaCycPathways = []
    keggPathways = []

    r = requests.get('http://' + dbDatabase + '.org/common/downloads/pathwayFiles//MetaCyc/')
    for row in r.text.split("\n"):
        if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
            metaCycPathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

    r = requests.get('http://' + dbDatabase + '.org/common/downloads/pathwayFiles//KEGG/')
    for row in r.text.split("\n"):
        if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
            keggPathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])
     
     if dbDatabase == 'tritrypdb':
        leishCycPathways = []
        trypanoCycPathways = []
        r = requests.get('http://' + dbDatabase + '.org/common/downloads/pathwayFiles//LeishCyc/')
        for row in r.text.split("\n"):
            if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
                leishCycPathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

        r = requests.get('http://' + dbDatabase + '.org/common/downloads/pathwayFiles//TrypanoCyc/')
        for row in r.text.split("\n"):
            if '<img src="/icons/unknown.gif" alt="[   ]"> <a href=' in row:
                trypanoCycPathways.append(row[len('<img src="/icons/unknown.gif" alt="[   ]"> <a href="'):].split('"')[0])

        dbDatabasePathways['LeishCyc'] = leishCycPathways
        dbDatabasePathways['TrypanoCyc'] = trypanoCycPathways

    dbDatabasePathways['MetaCyc'] = metaCycPathways
    dbDatabasePathways['KEGG'] = keggPathways

    return dbDatabasePathways

def linkPathwayIDrequestEuPathDB(dbDatabase, database, pathways):
    csvfile = open(temporaryDirectoryDatabase + "ecChebiToPathway_" + dbDatabase + database + ".tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(('pathway', 'ecChebis'))

    for pathway in pathways:
        if '.xgmml' in pathway:
            ecChebis = []
            r = requests.get('http://' + dbDatabase + '.org/common/downloads/pathwayFiles//' + database + '/' + pathway)
            for row in r.text.split("\n"):
                if "node label" in row:
                    ecChebis.append(row[len('  <node label="'):].split('"')[0])
            writer.writerow((pathway[:-len('.xgmml')], ecChebis))

    csvfile.close()

def main():
    dbDatabases = ['amoebadb', 'cryptodb', 'fungidb', 'giardiadb', 'microsporidiadb', 'piroplasmadb', 'plasmodb', 'toxodb', 'trichdb', 'tritrypdb']
    for dbDatabase in dbDatabases:
        dbDatabasePathways = requestEuPathDB(dbDatabase)
        linkPathwayIDrequestEuPathDB(dbDatabase, 'KEGG', dbDatabasePathways['KEGG'])
        linkPathwayIDrequestEuPathDB(dbDatabase, 'MetaCyc', dbDatabasePathways['MetaCyc'])
        if dbDatabase == 'tritrypdb':
            linkPathwayIDrequestEuPathDB(dbDatabase, 'LeishCyc', dbDatabasePathways['LeishCyc'])
            linkPathwayIDrequestEuPathDB(dbDatabase, 'TrypanoCyc', dbDatabasePathways['TrypanoCyc'])

main()
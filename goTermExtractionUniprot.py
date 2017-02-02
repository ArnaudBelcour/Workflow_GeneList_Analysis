#!/usr/bin/env python3

import csv
import math
import pandas as pa
import requests
import sys
from ast import literal_eval
from bs4 import BeautifulSoup
from functools import lru_cache
from SPARQLWrapper import SPARQLWrapper, JSON

inputDirectory = "inputFiles/"
temporaryDirectory = 'temporaryFiles/'
outputDirectory = 'outputFiles/'

def genomeFileExtraction():
    df = pa.read_csv(inputDirectory + "BrowserPlasmoUniprot.tsv", skiprows = 1, sep = "\t", header=None)
    df = df[[0, 1, 2, 3]]
    df.columns = [["GenePlasmo", "PosDepart", "PosFin", "Uniprot"]]

    #number = 299
    #listOfNumber = [0, 299]
    #for iterationOf708 in range(math.ceil(len(df["Uniprot"].tolist())/300)):
        #number = number + 300
        #listOfNumber.append(number)
    #del listOfNumber[-1]

    stringValuesProteinNames = []
    #stringValuesProteinName = ""

    #for mnemonicNameOfProtein in df["Uniprot"].tolist():
        #stringValuesProteinName = stringValuesProteinName + "'" + mnemonicNameOfProtein + "' "

    #for numberLimitationQuery in listOfNumber:
        #stringValuesProteinName = ""

        #if numberLimitationQuery != df["Uniprot"].tolist()[-1]:
            #currentList = df["Uniprot"].tolist()[numberLimitationQuery:numberLimitationQuery+354]
        #for mnemonicNameOfProtein in currentList:
            #stringValuesProteinName = "'" + mnemonicNameOfProtein + "'"
            #stringValuesProteinNames.append(stringValuesProteinName)

    for mnemonicNameOfProtein in df["Uniprot"].tolist():
        stringValuesProteinName = "'" + mnemonicNameOfProtein + "'"
        stringValuesProteinNames.append(stringValuesProteinName)

    return stringValuesProteinNames

def createApproximationGOGenome(stringValuesProteinNames):
    csvfile = open(temporaryDirectory + "test_genomeGO.tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(("mnemonicNameOfProtein", "uniprotID", 'GOs'))

    for stringValuesProteinName in stringValuesProteinNames:
        sparql = SPARQLWrapper("http://beta.sparql.uniprot.org/sparql")
        sparql.setQuery("""
        PREFIX up:<http://purl.uniprot.org/core/>

        SELECT ?protein ?go
        WHERE
        {
                ?protein a up:Protein .
                ?protein up:mnemonic ?proteinNamed .
                ?protein up:classifiedWith ?go .
                FILTER (regex(str(?go), "GO")) .
                VALUES ?proteinNamed {""" + stringValuesProteinName + """}
        }
        """)

        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

        goFounds = []
        for result in results["results"]["bindings"]:
            mnemonicNameOfProtein = stringValuesProteinName
            uniprotID = result["protein"]["value"][32:]
            goFounds.append(result["go"]["value"][31:])
        writer.writerow((mnemonicNameOfProtein, uniprotID, goFounds))

    csvfile.close()

    df = pa.read_csv(temporaryDirectory + 'test_genomeGO.tsv', '\t')

    for index, row in df.iterrows():
        row['GOs'] = unionGOAndTheirAncestors(literal_eval(row['GOs']))
    df.to_csv(temporaryDirectory + 'test_genomeGO_edit.tsv', '\t')

@lru_cache(maxsize = 12400)
def goTermAncestors(go):
    goAncestors = []
    sparql = SPARQLWrapper("http://localhost:3030/datanase/query")
    sparql.setQuery("""
    PREFIX rdfs:<http://www.w3.org/2000/01/rdf-schema#>
    PREFIX uniprot: <http://bio2rdf.org/uniprot:>
    PREFIX go: <http://purl.obolibrary.org/obo/GO_>
    PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
    PREFIX owl: <http://www.w3.org/2002/07/owl#>

    SELECT DISTINCT ?goAnc ?goAncLabel
    WHERE {
      go:""" + go[3:] + """ (rdfs:subClassOf|(rdfs:subClassOf/owl:someValuesFrom))*
    ?goAnc .
      OPTIONAL { ?goAnc rdfs:label ?goAncLabel .}
        FILTER ((str(?goAnc) != "b")) .
    }
    """)

    sparql.setReturnFormat(JSON)
    results = sparql.query().convert()

    for result in results["results"]["bindings"]:
        goAncestors.append(result["goAnc"]["value"][31:])

    return goAncestors

def unionGOAndTheirAncestors(gos):
    goAncestorsForGolists = []

    for go in gos:
        goAncestors = goTermAncestors(go)
        goAncestorsForGolists.append(goAncestors)

    golistForEntity = list(set().union(*goAncestorsForGolists))

    return golistForEntity

def requestGOTermForAProtein(uniprotID):
    r = requests.get('http://www.uniprot.org/uniprot/' + uniprotID + '.xml')
    soup = BeautifulSoup(r.text)
    l = soup.findAll({"type", "GO"} )

    td_tag_list = soup.findAll(
                    lambda tag:tag.name == "dbreference" and
                    tag["type"] == "GO")

    l_GOTerms = []

    for index in range(len(td_tag_list)):
        l_GOTerms.append(td_tag_list[index]["id"])

    return l_GOTerms

def inputPythonFormat(sentenceChoice, pythonVersion):
    if pythonVersion  < (3,0,0):
        choice = raw_input(sentenceChoice)
    if pythonVersion  > (3,0,0):
        choice = input(sentenceChoice)

    return choice

def main():
    #df = pa.read_csv("BrowserPlasmoUniprot.tsv", skiprows = 1, sep = "\t", header=None)
    #df = df[[0, 1, 2, 3]]
    #df.columns = [["GenePlasmo", "PosDepart", "PosFin", "Uniprot"]]

    #GOTerms = []

    #for uniprotID in df["Uniprot"].tolist():
        #GOTerms.extend(requestGOTermForAProtein(uniprotID))

    #df = pa.DataFrame(GOTerms)
    #df.columns = [["GOs"]]
    #df.to_csv(inputDirectory + "GOTermsPlasmoGenome.tsv", sep="\t", index = False, header = True, quoting = csv.QUOTE_NONE)
    pythonVersion = sys.version_info
    yesAnswers = ['y', 'yes']

    sentenceChoice = "Do you want to update the uniprotID from the genome? "
    yesOrNo = inputPythonFormat(sentenceChoice, pythonVersion)

    if yesOrNo.lower() in yesAnswers:
        stringValuesProteinNames = genomeFileExtraction()
        createApproximationGOGenome(stringValuesProteinNames)

main()

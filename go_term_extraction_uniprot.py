#!/usr/bin/env python3

import csv
import math
import pandas as pa
import requests
import six

from bs4 import BeautifulSoup
from functools import lru_cache
from SPARQLWrapper import SPARQLWrapper, JSON

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
output_directory = 'outputFiles/'

def genome_file_extraction():
    df = pa.read_csv(input_directory + "BrowserPlasmoUniprot.tsv", skiprows = 1, sep = "\t", header=None)
    df = df[[0, 1, 2, 3]]
    df.columns = [["GenePlasmo", "PosDepart", "PosFin", "Uniprot"]]

    #number = 299
    #listOfNumber = [0, 299]
    #for iterationOf708 in range(math.ceil(len(df["Uniprot"].tolist())/300)):
        #number = number + 300
        #listOfNumber.append(number)
    #del listOfNumber[-1]

    string_values_protein_names = []
    #string_values_protein_name = ""

    #for mnemonic_name_of_protein in df["Uniprot"].tolist():
        #string_values_protein_name = string_values_protein_name + "'" + mnemonic_name_of_protein + "' "

    #for numberLimitationQuery in listOfNumber:
        #string_values_protein_name = ""

        #if numberLimitationQuery != df["Uniprot"].tolist()[-1]:
            #currentList = df["Uniprot"].tolist()[numberLimitationQuery:numberLimitationQuery+354]
        #for mnemonic_name_of_protein in currentList:
            #string_values_protein_name = "'" + mnemonic_name_of_protein + "'"
            #string_values_protein_names.append(string_values_protein_name)

    for mnemonic_name_of_protein in df["Uniprot"].tolist():
        string_values_protein_name = "'" + mnemonic_name_of_protein + "'"
        string_values_protein_names.append(string_values_protein_name)

    return string_values_protein_names

def create_approximation_go_genome(string_values_protein_names):
    csvfile = open(temporary_directory + "test_genomeGO.tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(("mnemonic_name_of_protein", "uniprot_id", 'GOs'))

    for string_values_protein_name in string_values_protein_names:
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
                VALUES ?proteinNamed {""" + string_values_protein_name + """}
        }
        """)

        sparql.setReturnFormat(JSON)
        results = sparql.query().convert()

        go_founds = []
        for result in results["results"]["bindings"]:
            mnemonic_name_of_protein = string_values_protein_name
            uniprot_id = result["protein"]["value"][32:]
            go_founds.append(result["go"]["value"][31:])
        writer.writerow((mnemonic_name_of_protein, uniprot_id, go_founds))

    csvfile.close()

    df = pa.read_csv(temporary_directory + 'test_genomeGO.tsv', '\t')

    for index, row in df.iterrows():
        row['GOs'] = union_go_and_their_ancestor(row['GOs'].split(","))
    df.to_csv(temporary_directory + 'test_genomeGO_edit.tsv', '\t')

def request_go_term_for_a_protein(uniprot_id):
    '''
        Request Uniprot to extract GO terms associated with a Uniprot ID.
    '''
    r = requests.get('http://www.uniprot.org/uniprot/' + uniprot_id + '.xml')
    soup = BeautifulSoup(r.text)
    l = soup.findAll({"type", "GO"} )

    td_tag_list = soup.findAll(
                    lambda tag:tag.name == "dbreference" and
                    tag["type"] == "GO")

    go_terms = []

    for index in range(len(td_tag_list)):
        go_terms.append(td_tag_list[index]["id"])

    return go_terms

def main():
    #df = pa.read_csv("BrowserPlasmoUniprot.tsv", skiprows = 1, sep = "\t", header=None)
    #df = df[[0, 1, 2, 3]]
    #df.columns = [["GenePlasmo", "PosDepart", "PosFin", "Uniprot"]]

    #GOTerms = []

    #for uniprot_id in df["Uniprot"].tolist():
        #GOTerms.extend(requestGOTermForAProtein(uniprot_id))

    #df = pa.DataFrame(GOTerms)
    #df.columns = [["GOs"]]
    #df.to_csv(input_directory + "GOTermsPlasmoGenome.tsv", sep="\t", index = False, header = True, quoting = csv.QUOTE_NONE)
    yes_answers = ['y', 'yes']

    yes_or_no = input("Do you want to update the uniprot_id from the genome? ")

    if yes_or_no.lower() in yes_answers:
        string_values_protein_names = genome_file_extraction()
        create_approximation_go_genome(string_values_protein_names)

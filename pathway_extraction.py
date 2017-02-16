#!/usr/bin/env python3

import math
import os
import pandas as pa
import re
import subprocess
from ast import literal_eval

temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'

def ipr_extraction(df_genome):
    ipr_expression = r"IPR[\d]{6}[^0-9]"

    df_genome = df_genome[['Gene_Name', 'GOs', 'EnzymeCodes', 'InterProScan']]
    df_genome = df_genome.set_index("Gene_Name")
    df_genome['InterProScan'] = df_genome['InterProScan'].str.split("; ").apply(
        lambda x: [y[:len('IPRXXXXXX')] 
                   for y in x 
                   if re.match(ipr_expression, y)])

def ec_extraction(df_genome):
    df_genome['EnzymeCodes'] = df_genome['EnzymeCodes'].str.lower().str.split(";")

    ecs_requests = []
    for ecs in df_genome['EnzymeCodes']:
        if type(ecs) != float:
            for ec in ecs:
                if ec not in ecs_requests:
                    ecs_requests.append(ec)

    return ecs_requests

def r_keggrest_ec(ecs_requests):
    command = 'Rscript'
    path_script = 'pathway-extraction/enzyme_to_pathway.R'
    cmd = [command, path_script] + ecs_requests
    try:
        subprocess.check_output(cmd, universal_newlines=True,stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

        temporary_directory_database = "temporaryFiles/databases/"

def translation_go_chebi(selected_gos):
    chebis = []
    for selected_go in literal_eval(selected_gos):
        if selected_go in go_to_chebi.index:
            if type(go_to_chebi.loc[selected_go]['ChEBI']) == str :
                chebis.append(go_to_chebi.loc[selected_go]['ChEBI'])
            if type(go_to_chebi.loc[selected_go]['ChEBI']) == list :
                chebis.extend(go_to_chebi.loc[selected_go]['ChEBI'].tolist())

    return set(chebis)

def main():
    if os.path.exists(temporary_directory_database[:-1]) == False :
        os.makedirs(temporary_directory_database)
    name_reference_file = 'Annotation_blast2go_PROT_eH'
    df_genome = pa.read_csv(temporary_directory + name_reference_file + 'GOsTranslatedAndFixed.tsv', sep = "\t")
    ecs_requests = ec_extraction(df_genome)
    r_keggrest_ec(ecs_requests)

    go_to_chebi = pa.read_csv(temporary_directory_database + "go_chebi_mapping.tsv", sep = '\t')
    go_to_chebi = go_to_chebi.set_index("GOs")
    df_genome['ChEBI'] = df_genome['GOs'].apply(translation_go_chebi)
main()

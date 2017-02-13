#!/usr/bin/env python3

import math
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
    path_script = 'pathway-extraction/enzymeToPathway.R'
    cmd = [command, path_script] + ecs_requests
    try:
        subprocess.check_output(cmd, universal_newlines=True,stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

def main():
    if os.path.exists(temporary_directory_database[:-1]) == False :
        os.makedirs(temporary_directory_database)
    name_reference_file = 'Annotation_blast2go_PROT_eHGOs'
    df_genome = pa.read_csv(temporary_directory + name_reference_file + 'GOsTranslatedAndFixed.tsv', sep = "\t")
    ecs_requests = ec_extraction(df_genome)
    r_keggrest_ec(ecs_requests)

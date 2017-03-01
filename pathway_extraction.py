#!/usr/bin/env python3

import os
import pandas as pa
import re
import subprocess

temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'

def ipr_extraction(df_genome):
    ipr_expression = r"IPR[\d]{6}[^0-9]"

    df_genome = df_genome[['Gene_Name', 'GOs', 'EnzymeCodes', 'InterProScan']]
    df_genome.set_index("Gene_Name", inplace = True)
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
    '''
        Uses the R script to extract pathway using Enzyme Code in KEGG.
    '''
    command = 'Rscript'
    path_script = 'pathway_extraction/keggrest_pathway_extraction.R'
    data_name = ["enzyme"]
    cmd = [command, path_script] + ecs_requests + data_name
    try:
        subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

def translation_go_data(selected_datas, df_mapping, data_column):
    datas = []
    for selected_data in selected_datas.split(","):
        if selected_data in df_mapping.index:
            if type(df_mapping.loc[selected_data][data_column]) == str :
                datas.append(df_mapping.loc[selected_data][data_column])
            if type(df_mapping.loc[selected_data][data_column]) == list :
                datas.extend(df_mapping.loc[selected_data][data_column].tolist())

    return set(datas)

def mapping_data(file_name, df_genome):
    df_mapping = pa.read_csv(temporary_directory_database + file_name, sep = '\t')
    df_mapping.set_index("GOs", inplace = True)
    data_column = df_mapping.columns[0]
    df_genome[data_column] = df_genome['GOs'].apply(translation_go_data, args = (df_mapping, data_column))

    return df_genome

def main():
    if os.path.exists(temporary_directory_database[:-1]) == False :
        os.makedirs(temporary_directory_database)
    name_reference_file = 'Annotation_blast2go_PROT_eH'
    df_genome = pa.read_csv(temporary_directory + name_reference_file + 'GOsTranslatedAndFixed.tsv', sep = "\t")
    ecs_requests = ec_extraction(df_genome)
    r_keggrest_ec(ecs_requests)

    for file_name in os.listdir(temporary_directory_database):
        if "mapping" in file_name:
            df_genome = mapping_data(file_name, df_genome)

main()

#!/usr/bin/env python3

import numpy as np
import os
import pandas as pa

from ast import literal_eval

input_directory = "inputFiles/"
temporary_directory = "temporaryFiles/"
temporary_directory_database = "temporaryFiles/databases/"

def translation_data(selected_datas, df_mapping, data_column, specification):
    datas = []
    for selected_data in selected_datas.split(","):
        if selected_data in df_mapping.index:
            data_retrieved = df_mapping.loc[selected_data][data_column]
            if type(data_retrieved) == list :
                datas.extend(data_retrieved.tolist())
            elif type(data_retrieved) == str :
                datas.append(str(data_retrieved))
    datas = list(set(datas))
    if specification == 'pathway':
        return datas
    elif specification == 'initial':
        datas_string = ",".join(datas)
        return datas_string

def translation_interpro_data(selected_interpros, df_mapping, data_column, database):
    datas = []

    for selected_interpro in selected_interpros.split(','):
        if selected_interpro in df_mapping.index:
            if type(df_mapping['Database'].loc[selected_interpro]) is list:
                if all(df_mapping['Database'].loc[selected_interpro].tolist()) == database:
                    if type(df_mapping.loc[selected_interpro][data_column]) not in [str, float]:
                        datas.extend(df_mapping.loc[selected_interpro][data_column].tolist())
                    elif type(df_mapping.loc[selected_interpro][data_column]) is str or type(df_mapping.loc[selected_interpro][data_column]) is float:
                        datas.append(str(df_mapping.loc[selected_interpro][data_column]))
            if type(df_mapping['Database'].loc[selected_interpro]) is str:
                if df_mapping['Database'].loc[selected_interpro] == database:
                    if type(df_mapping.loc[selected_interpro][data_column]) not in [str, float]:
                        datas.extend(df_mapping.loc[selected_interpro][data_column].tolist())
                    elif type(df_mapping.loc[selected_interpro][data_column]) is str or type(df_mapping.loc[selected_interpro][data_column]) is float:
                        datas.append(str(df_mapping.loc[selected_interpro][data_column]))

    datas = list(set(datas))

    return datas

def translation_gene_pathway(selected_gene, df_mapping, data_column):
    datas = []

    if selected_gene in df_mapping.index:
        if type(df_mapping.loc[selected_gene][data_column]) is str:
            datas = literal_eval(df_mapping.loc[selected_gene][data_column])

    datas = list(set(datas))

    return datas

def intialize_column(df, column_name):
    df[column_name] = ''
    df[column_name] = df[column_name].apply(list)

    return df

def drop_duplicates(datas):
    datas_with_unique = []
    for data in datas:
        if data not in datas_with_unique:
            datas_with_unique.append(data)
    return datas_with_unique

def list_to_string(datas):
    return ','.join(datas)

def main(file_name_temporary):
    df_genome = pa.read_csv(temporary_directory + file_name_temporary, sep = "\t")
    df_genome.replace(np.nan, '', regex=True, inplace=True)

    for file_name in os.listdir(temporary_directory_database):
        if "ecChebiToPathway" in file_name:
            df_eupathdb = pa.read_csv(temporary_directory_database + file_name, sep = '\t')
            df_eupathdb['ecChebis'] = df_eupathdb['ecChebis']
            df_eupathdb = df_eupathdb.set_index('ecChebis')
            df_genome['pathway_' + file_name] = (df_genome['EnzymeCodes'].apply(translation_data, args = (df_eupathdb, 'pathway', 'pathway'))\
                                                 + df_genome['ChEBI'].apply(translation_data, args = (df_eupathdb, 'pathway', 'pathway')))

    column_names = df_genome.columns.tolist()

    kegg_eupathdb_columns = []
    metacyc_eupathdb_columns = []
    leishcyc_eupathdb_columns = []
    trypanocyc_eupathdb_columns = []

    for column_name in column_names:
        if 'ecChebiToPathway' in column_name and 'KEGG' in column_name:
            kegg_eupathdb_columns.append(column_name)
        if 'ecChebiToPathway' in column_name and 'MetaCyc' in column_name:
            metacyc_eupathdb_columns.append(column_name)
        if 'ecChebiToPathway' in column_name and 'LeishCyc' in column_name:
            leishcyc_eupathdb_columns.append(column_name)
        if 'ecChebiToPathway' in column_name and 'TrypanoCyc' in column_name:
            trypanocyc_eupathdb_columns.append(column_name)

    for database in ['kegg', 'metacyc', 'leishcyc', 'trypanocyc']:
        df_genome = intialize_column(df_genome, 'pathway_eupathdb_' + database)

    for column in kegg_eupathdb_columns:
        df_genome['pathway_eupathdb_kegg'] += df_genome[column]
    for column in metacyc_eupathdb_columns:
        df_genome['pathway_eupathdb_metacyc'] += df_genome[column]
    for column in leishcyc_eupathdb_columns:
        df_genome['pathway_eupathdb_leishcyc'] += df_genome[column]
    for column in trypanocyc_eupathdb_columns:
        df_genome['pathway_eupathdb_trypanocyc'] += df_genome[column]

    for column_name in column_names:
        if 'ecChebiToPathway' in column_name:
            df_genome = df_genome.drop(column_name, 1)

    df_genome['pathway_eupathdb_kegg'] = df_genome['pathway_eupathdb_kegg'].apply(drop_duplicates)
    df_genome['pathway_eupathdb_metacyc'] = df_genome['pathway_eupathdb_metacyc'].apply(drop_duplicates)
    df_genome['pathway_eupathdb_leishcyc'] = df_genome['pathway_eupathdb_leishcyc'].apply(drop_duplicates)
    df_genome['pathway_eupathdb_trypanocyc'] = df_genome['pathway_eupathdb_trypanocyc'].apply(drop_duplicates)

    df_genome['pathway_eupathdb_kegg'] = df_genome['pathway_eupathdb_kegg'].apply(list_to_string)
    df_genome['pathway_eupathdb_metacyc'] = df_genome['pathway_eupathdb_metacyc'].apply(list_to_string)
    df_genome['pathway_eupathdb_leishcyc'] = df_genome['pathway_eupathdb_leishcyc'].apply(list_to_string)
    df_genome['pathway_eupathdb_trypanocyc'] = df_genome['pathway_eupathdb_trypanocyc'].apply(list_to_string)

    for file_name in os.listdir(temporary_directory_database):
        if "pathway_reactome" in file_name:
            df_reactome = pa.read_csv(temporary_directory_database + file_name, sep = "\t")
            if 'GO' in file_name:
                df_reactome = df_reactome.set_index('GO')
                df_genome[file_name] = df_genome['GOs'].apply(translation_data, args = (df_reactome, 'Id', 'pathway'))
            if 'Interpro' in file_name:
                df_reactome = df_reactome.set_index('Interpro')
                df_genome[file_name] = df_genome['InterProScan'].apply(translation_data, args = (df_reactome, 'Id', 'pathway'))
            if 'EC' in file_name:
                df_reactome['EC'] = df_reactome['EC'].str.replace("ec:", "")
                df_reactome['EC'] = df_reactome['EC'].str.replace("_", ":")
                df_reactome = df_reactome.set_index('EC')
                df_genome[file_name] = df_genome['EnzymeCodes'].apply(translation_data, args = (df_reactome, 'Id', 'pathway'))
            if 'CHEBI' in file_name:
                df_reactome['CHEBI'] = df_reactome['CHEBI']
                df_reactome = df_reactome.set_index('CHEBI')
                df_genome[file_name] = df_genome['ChEBI'].apply(translation_data, args = (df_reactome, 'Id', 'pathway'))

    column_names = df_genome.columns.tolist()
    pathway_reactomes = []
    for column_name in column_names:
        if 'pathway_reactome' in column_name:
            pathway_reactomes.append(column_name)

    df_genome = intialize_column(df_genome, 'pathway_reactome')

    for column in pathway_reactomes:
        df_genome['pathway_reactome'] += df_genome[column]

    for column_name in column_names:
        if 'pathway_reactome_' in column_name:
            df_genome = df_genome.drop(column_name, 1)

    df_genome['pathway_reactome'] = df_genome['pathway_reactome'].apply(drop_duplicates)
    df_genome['pathway_reactome'] = df_genome['pathway_reactome'].apply(list_to_string)

    df_keggrest = pa.read_csv(temporary_directory_database + 'enzyme_pathway_kegg.tsv', sep='\t')
    df_keggrest['ecCode'] = df_keggrest['ecCode'].str.replace("ec:", "")
    df_keggrest = df_keggrest.set_index('ecCode')

    df_genome['pathway_keggrest'] = df_genome['EnzymeCodes'].apply(translation_data, args = (df_keggrest, 'PathwayID', 'pathway'))

    df_genome['pathway_keggrest'] = df_genome['pathway_keggrest'].apply(list_to_string)

    df_interpro = pa.read_csv(temporary_directory_database + 'interpro_pathway.tsv', sep='\t')
    df_interpro = df_interpro.set_index('Interpros')

    databases = ['KEGG', 'REACTOME', 'METACYC']

    for database in databases:
        df_genome['pathway_interpro_' + database] = df_genome['InterProScan'].apply(translation_interpro_data, args = (df_interpro, 'Pathway_id', database))
        df_genome['pathway_interpro_' + database] = df_genome['pathway_interpro_' + database].apply(list_to_string)


    df_ghost_koala = pa.read_csv(temporary_directory_database + 'gene_with_kegg_pathway.tsv', sep='\t')
    df_ghost_koala = df_ghost_koala.set_index('Gene')

    df_genome['pathway_ghost_koala'] = df_genome['Gene_Name'].apply(translation_gene_pathway, args = (df_ghost_koala, 'kegg_pathway'))
    df_genome['pathway_ghost_koala'] = df_genome['pathway_ghost_koala'].apply(drop_duplicates)
    df_genome['pathway_ghost_koala'] = df_genome['pathway_ghost_koala'].apply(list_to_string)

    df_genome.to_csv(temporary_directory + "result_pathway_extraction.tsv", sep="\t", index=False)


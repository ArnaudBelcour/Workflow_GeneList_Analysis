#!/usr/bin/env python3

import csv
import math
import numpy as np
import os
import pandas as pa
import pronto
import re
import requests
import shutil
import six

from tqdm import *

import ancestor_go_extraction
import mapping_pathway_data
import pathway_extractor
import pathway_extraction.uniprot_retrieval_data as uniprot_retrieval_data

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'

class FileManagement():

    def __init__(self, name_of_the_file):
        self._file_name, self._file_extension = os.path.splitext(name_of_the_file)

    @property
    def file_name(self):
        return self._file_name

    @file_name.setter
    def file_name(self, name):
        self._file_name = name

    @property
    def file_extension(self):
        return self._file_extension

    @file_extension.setter
    def file_extension(self, extension):
        self._file_extension = extension

    def go_label_number_dictionary_creation(self, specification='normal'):
        '''
            Check the internet connection to use a http request or a file.
        '''
        try:
            requests.get('http://www.google.com/')
            return self.go_label_number_dictionary_creation_from_http(specification)

        except requests.ConnectionError:
            return self.go_label_number_dictionary_creation_from_file(specification)

    def go_label_number_dictionary_creation_from_http(self, specification='normal'):
        '''
            Create a dictionary containing GO labels (as key) associated with their GO numbers (as value), if specification is 'normal'. Use to translate GO labels into GO numbers.
            Or GO numbers (as key) associated with their GO labels (as value), if specification is 'inverse'. Use to translate GO numbers into GO labels.
            Create also a dictionary containing the synonym of the GO terms.
        '''
        def dict_to_file(dictionary, file_name, specification='normal'):
            df = pa.DataFrame.from_dict(dictionary, orient='index')
            df.reset_index(inplace=True)
            if specification == 'inverse':
                df.columns = [['GOnumber', 'GOlabel']]
            elif specification == 'normal':
                df.columns = [['GOlabel', 'GOnumber']]
            df.to_csv(temporary_directory + file_name + '.tsv', sep='\t', index=False)

        d_go_label_to_number = {}

        go_ontology = pronto.Ontology('http://purl.obolibrary.org/obo/go/go-basic.obo')

        if specification == "inverse":
            for go_term in go_ontology:
                d_go_label_to_number[go_term.id] = go_term.name
            dict_to_file(d_go_label_to_number, 'go_number_label', specification)

            return d_go_label_to_number

        else:
            for go_term in go_ontology:
                d_go_label_to_number[go_term.name] = go_term.id

        d_go_label_with_synonym = {}

        for term in go_ontology:
            if getattr(term, 'synonyms', 'default value'):
                for synonym in term.synonyms:
                    d_go_label_with_synonym[str(synonym).split('"')[1]] = term.id

        dict_to_file(d_go_label_to_number, 'go_number_label')
        dict_to_file(d_go_label_with_synonym, 'go_number_label_synonym')

        return d_go_label_to_number, d_go_label_with_synonym

    def go_label_number_dictionary_creation_from_file(self, specification='normal'):
        if specification == "inverse":
            df = pa.read_csv(temporary_directory + 'go_number_label.tsv', sep='\t')
            df.set_index('GOnumber', inplace=True)
            d_go_label_to_number = df.to_dict('dict')['GOlabel']

            return d_go_label_to_number

        else:
            df = pa.read_csv(temporary_directory + 'go_number_label.tsv', sep='\t')
            df.set_index('GOlabel', inplace=True)
            d_go_label_to_number = df.to_dict('dict')['GOnumber']

        df = pa.read_csv(temporary_directory + 'go_number_label_synonym.tsv', sep='\t')
        df.set_index('GOlabel', inplace=True)
        d_go_label_with_synonym = df.to_dict('dict')['GOnumber']

        return d_go_label_to_number, d_go_label_with_synonym

    def find_column_of_interest(self, df):
        '''
            Detect columns containing GO number, GO label, EC number and Interpro ID.
            To do this, regular expression are used, for each types of data.
            The occurrence of each regular expression is counted.
            Then the column containing the maximum of occurrence for a type of data is associated with it by returning it's name.
            Between GO labels and GO numbers, GO numbers are preferred.
        '''
        columns = df.columns.tolist()

        go_label_expression = r"[FPC]{1}:[\w]*"
        go_number_expression = r"[FPC]{1}:GO[:_][\d]{7}"
        ec_expression = r"[Ee][Cc]:[\d]{1}[\.]{1}[\d]{,2}[\.]{,1}[\d]{,2}[\.]{,1}[\d]{,2}"
        ipr_expression = r"IPR[\d]{6}"
        ko_kegg_expression = r"K[\d]{5}"

        go_label_columns = {}
        go_number_columns = {}
        ec_columns = {}
        ipr_columns = {}
        ko_keggs = {}

        for column in columns:
            for column_values in df[column]:
                if type(column_values) is not str:
                    continue
                if re.match(go_number_expression, column_values):
                    if column in go_number_columns:
                        go_number_columns[column] += 1
                    else:
                        go_number_columns[column] = 1
                elif re.match(go_label_expression, column_values):
                    if column in go_label_columns:
                        go_label_columns[column] += 1
                    else:
                        go_label_columns[column] = 1
                elif re.match(ec_expression, column_values):
                    if column in ec_columns:
                        ec_columns[column] += 1
                    else:
                        ec_columns[column] = 1
                elif re.match(ipr_expression, column_values):
                    if column in ipr_columns:
                        ipr_columns[column] += 1
                    else:
                        ipr_columns[column] = 1
                elif re.match(ko_kegg_expression, column_values):
                    if column in ko_keggs:
                        ko_keggs[column] += 1
                    else:
                        ko_keggs[column] = 1

        if go_number_columns:
            go_number_column = max(go_number_columns, key=go_number_columns.get)
            go_column = go_number_column
            go_label_columns.pop(go_number_column, None)
        if ec_columns != []:
            ec_column = max(ec_columns, key=ec_columns.get)
        else:
            ec_column = np.nan
        if ipr_columns != []:
            ipr_column = max(ipr_columns, key=ipr_columns.get)
        else:
            ipr_column = np.nan
        if ko_keggs and ko_keggs != []:
            ko_kegg = max(ko_keggs, key=ko_keggs.get)
        else:
            ko_kegg = np.nan

        if not go_number_columns:
            go_label_column = max(go_label_columns, key=go_label_columns.get)
            go_column = go_label_column

        return go_column, ec_column, ipr_column

    def preprocessing_file(self, type_file):

        def drop_duplicates(datas):
            datas_with_unique = []
            for data in datas:
                if data not in datas_with_unique:
                    datas_with_unique.append(data)

            return datas_with_unique

        name_input_file = self.file_name
        extension_input_file = self.file_extension

        if extension_input_file == '.xls':
            results_dataframe = pa.read_excel(input_directory + name_input_file + extension_input_file, sep=None, na_values="")
        else:
                #Check if the file contains more than one column (using sep=None, checking some delimiters) if not a TypeError occurres (there is only one column).
            try:
                results_dataframe = pa.read_csv(input_directory + name_input_file + extension_input_file, sep=None, engine="python", na_values="")
            except TypeError:
                results_dataframe = pa.read_csv(input_directory + name_input_file + extension_input_file, header=None, engine="python")

        yes_or_no = input("Is the first columns of your file, the column containing gene name? ")

        yes_answers = ['yes', 'y', 'oui', 'o']
        if yes_or_no.lower() in yes_answers :
            name_gene_column = results_dataframe.columns[0]
        else :
            name_gene_column = input("Write the name of the column containing the gene names : ")

        if type_file == "gene_list":
            results_dataframe = results_dataframe[[name_gene_column]]
            results_dataframe.columns = [['Gene_Name']]
            results_dataframe.to_csv(temporary_directory + name_input_file + "GOsTranslatedAndFixed.tsv", "\t", index=False, header=True, quoting=csv.QUOTE_NONE)

            return name_input_file + "GOsTranslatedAndFixed.tsv"

        go_column, ec_column, ipr_column = self.find_column_of_interest(results_dataframe)

        uniprot_retrieval_y_n = input("Do you want to try to retrieve data from blast results? ").lower()
        if uniprot_retrieval_y_n in yes_answers :
            column_name = input("What is the name of the column of the blast results? ")
            results_dataframe = results_dataframe[[name_gene_column, go_column, ec_column, ipr_column, column_name]]
            results_dataframe.columns = [['Gene_Name', 'GOs', 'EnzymeCodes', 'InterProScan', 'Blast']]
        else:
            results_dataframe = results_dataframe[[name_gene_column, go_column, ec_column, ipr_column]]
            results_dataframe.columns = [['Gene_Name', 'GOs', 'EnzymeCodes', 'InterProScan']]

        results_dataframe['GOs'] = results_dataframe['GOs'].str.replace("C:", "")
        results_dataframe['GOs'] = results_dataframe['GOs'].str.replace("P:", "")
        results_dataframe['GOs'] = results_dataframe['GOs'].str.replace("F:", "")
        results_dataframe['GOs'] = results_dataframe['GOs'].str.split("; ")
        results_dataframe['GOs'] = results_dataframe['GOs'].str.join(",")

        results_dataframe.replace(np.nan, '', regex=True, inplace=True)
        list_to_string = lambda datas: ','.join(datas)
        to_list = lambda x: x.split("; ")

        results_dataframe['InterProScan'] = results_dataframe['InterProScan'].apply(to_list)

        ipr_expression = r"IPR[\d]{6}"
        ipr_selection = lambda interpros: [interpro[:9]
                                           for interpro in interpros
                                           if re.match(ipr_expression, interpro[:9])]

        results_dataframe['InterProScan'] = results_dataframe['InterProScan'].apply(ipr_selection)
        results_dataframe['InterProScan'] = results_dataframe['InterProScan'].apply(drop_duplicates)
        results_dataframe['InterProScan'] = results_dataframe['InterProScan'].apply(list_to_string)

        results_dataframe['EnzymeCodes'] = results_dataframe['EnzymeCodes'].apply(to_list)

        ec_modification = lambda ecs: [ec.lower()
                                     for ec in ecs]

        results_dataframe['EnzymeCodes'] = results_dataframe['EnzymeCodes'].apply(ec_modification)
        results_dataframe['EnzymeCodes'] = results_dataframe['EnzymeCodes'].apply(drop_duplicates)
        results_dataframe['EnzymeCodes'] = results_dataframe['EnzymeCodes'].apply(list_to_string)

        if uniprot_retrieval_y_n in yes_answers :
            results_dataframe = uniprot_retrieval_data.extract_information_from_uniprot(results_dataframe)

        results_dataframe.to_csv(temporary_directory + name_input_file + "GOsTranslatedAndFixed.tsv", "\t", index=False, quoting=csv.QUOTE_NONE)

        return name_input_file + "GOsTranslatedAndFixed.tsv"

class FileManagementGeneGO(FileManagement):

    def __init__(self, name_of_the_file, type_of_the_file, column_name):
        FileManagement.__init__(self, name_of_the_file)
        self._analyzed_object = column_name
        self._type_file = type_of_the_file

    @property
    def analyzed_object_name(self):
        return self._analyzed_object

    @analyzed_object_name.setter
    def analyzed_object_name(self, name):
        self._analyzed_object = name

    @property
    def type_file(self):
        return self._type_file

    @type_file.setter
    def type_file(self, type_name):
        self._type_file = type_name

    def go_ancestors_list_of_interest(self, column_analyzed_object, file_name_temporary):

        print("GO ancestors retrieval")
        df = pa.read_csv(temporary_directory + file_name_temporary, sep='\t')
        df.replace(np.nan, '', regex=True, inplace=True)

        for index, row in tqdm(df.iterrows(), total=len(df.index)):
            df[column_analyzed_object].loc[index] = ancestor_go_extraction.union_go_and_their_ancestor(row[column_analyzed_object].split(","))

        df.to_csv(temporary_directory + file_name_temporary, sep='\t', index=False)

    def counting_genome(self, file_name_temporary, column_name, column_analyzed_object):
        analyzed_objects = []
        df = pa.read_csv(temporary_directory + file_name_temporary, sep="\t")
        df.replace(np.nan, '', regex=True, inplace=True)
        df.set_index('Gene_Name', inplace=True)

        for index, row in df.iterrows():
            if row['GOs'] == '':
                df.drop(index, inplace=True)

        for index, row in df.iterrows():
            for analyzed_object in row[column_analyzed_object].split(","):
                analyzed_objects.append(analyzed_object)

        counts_df_genome = pa.DataFrame(analyzed_objects)
        counts_df_genome.columns = [column_analyzed_object]
        counts_df_genome = counts_df_genome.groupby(column_analyzed_object).size().rename(column_name)
        counts_df_genome = counts_df_genome.to_frame()
        counts_df_genome.reset_index(inplace=True)
        counts_df_genome.to_csv(temporary_directory + "counting_objects_in_genome.tsv", "\t", index=False, header=True, quoting=csv.QUOTE_NONE)

        number_of_gene = len(df.index.unique())

        return "counting_objects_in_genome", number_of_gene

class FileManagementGeneGOsGenome(FileManagementGeneGO):

    def __init__(self, name_of_the_file, type_of_the_file, column_name):
        FileManagementGeneGO.__init__(self, name_of_the_file, type_of_the_file, column_name)

    def genome_file_processing(self):
        string_to_boolean = lambda value: value.lower() in ("yes", "true", "y", "t", "1")

        analyzed_object_name = self.analyzed_object_name
        name_of_the_file = self.file_name
        extension_of_the_file  = self.file_extension
        type_of_the_file = self.type_file

        already_analyzed_file_yes_no = string_to_boolean(input("Does this file have already been analyzed? "))
        if already_analyzed_file_yes_no == True:
            shutil.copy(input_directory + name_of_the_file + extension_of_the_file, temporary_directory)
            file_name_temporary = name_of_the_file + extension_of_the_file

        elif already_analyzed_file_yes_no == False:
            file_name_temporary = self.preprocessing_file(type_of_the_file)

            request_base_yes_no = string_to_boolean(input("Do you want to update your databases? "))
            if request_base_yes_no == True:
                session = requests.Session()
                pathway_extractor.data_retrieval_from_GO(file_name_temporary)
                pathway_extractor.main(file_name_temporary, session)

            self.go_ancestors_list_of_interest(analyzed_object_name, file_name_temporary)
            mapping_pathway_data.main(file_name_temporary)

        counting_object_file, number_of_gene_genome = self.counting_genome(file_name_temporary, 'CountsReference', analyzed_object_name)

        return file_name_temporary, counting_object_file, number_of_gene_genome

class FileManagementGeneGOGenome(FileManagementGeneGO):

    def __init__(self, name_of_the_file, type_of_the_file, column_name):
        FileManagementGeneGO.__init__(self, name_of_the_file, type_of_the_file, column_name)

    def genome_file_processing(self):
        analyzed_object_name = self.analyzed_object_name
        name_of_the_file = self.file_name
        extension_of_the_file  = self.file_extension

        df = pa.read_csv(input_directory + name_of_the_file + extension_of_the_file, sep="\t", header=None)

        df.columns = [['Gene_Name', 'GOs']]
        df = df.groupby('Gene_Name')['GOs'].apply(';'.join)
        df = df.to_frame()

        df['GOs'] = df['GOs'].str.replace("C:", "")
        df['GOs'] = df['GOs'].str.replace("P:", "")
        df['GOs'] = df['GOs'].str.replace("F:", "")
        df['GOs'] = df['GOs'].str.split(";")
        df['GOs'] = df['GOs'].str.join(",")

        df.replace(np.nan, '', regex=True, inplace=True)

        file_name_temporary = name_of_the_file + "GOsTranslatedAndFixed.tsv"

        df.to_csv(temporary_directory + file_name_temporary, "\t", index=True, quoting=csv.QUOTE_NONE)

        self.go_ancestors_list_of_interest(analyzed_object_name, file_name_temporary)

        counting_object_file, number_of_gene = self.counting_genome(file_name_temporary, 'CountsReference', analyzed_object_name)

        return file_name_temporary, counting_object_file, number_of_gene

class FileManagementGeneInterest(FileManagementGeneGO):

    def __init__(self, name_of_the_file, type_of_the_file, column_name, file_name_genome):
        FileManagementGeneGO.__init__(self, name_of_the_file, type_of_the_file, column_name)
        self._file_genome_reference_name = file_name_genome

    @property
    def genome_file_reference_name(self):
        return self._file_genome_reference_name

    @genome_file_reference_name.setter
    def genome_file_reference_name(self, file_name):
        self._file_genome_reference_name = file_name

    def interest_file_processing(self):
        string_to_boolean = lambda value: value.lower() in ("yes", "true", "y", "t", "1")

        analyzed_object_name = self.analyzed_object_name
        name_of_the_file = self.file_name
        extension_of_the_file  = self.file_extension
        type_of_the_file = self.type_file

        file_name_temporary = self.preprocessing_file(type_of_the_file)
        counting_object_file, number_of_gene_list = self.counting_gene_list(file_name_temporary, 'Counts', analyzed_object_name)

        return counting_object_file, number_of_gene_list

    def counting_gene_list(self, file_name_temporary, column_name, column_analyzed_object):
        analyzed_objects = []
        df = pa.read_csv(temporary_directory + file_name_temporary, sep="\t")
        df = df[['Gene_Name']]
        df.set_index("Gene_Name", inplace=True)

        df_genome = pa.read_csv(temporary_directory + self.genome_file_reference_name, sep="\t")
        df_genome = df_genome[['Gene_Name', 'GOs']]
        df_genome.set_index("Gene_Name", inplace=True)

        df_joined = df.join(df_genome)
        df_joined = df_joined[pa.notnull(df_joined[column_analyzed_object])]
        df_joined.reset_index(inplace=True)
        df_joined.to_csv(temporary_directory + file_name_temporary, "\t", index=False, header=True, quoting=csv.QUOTE_NONE)
        df_joined.set_index("Gene_Name", inplace=True)

        for index, row in df_joined.iterrows():
            if row['GOs'] == '':
                df_joined.drop(index, inplace=True)

        for index, row in df_joined.iterrows():
            for analyzed_object in row[column_analyzed_object].split(","):
                analyzed_objects.append(analyzed_object)

        counts_df = pa.DataFrame(analyzed_objects)
        counts_df.columns = [column_analyzed_object]
        counts_df = counts_df.groupby(column_analyzed_object).size().rename(column_name)
        counts_df = counts_df.to_frame()

        number_of_gene = len(df.index.unique())

        counts_df.reset_index(inplace=True)
        counts_df.to_csv(temporary_directory + "counting_objects_in_interest.tsv", "\t", index=False, header=True, quoting=csv.QUOTE_NONE)

        return "counting_objects_in_interest", number_of_gene

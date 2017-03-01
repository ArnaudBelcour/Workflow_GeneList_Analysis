#!/usr/bin/env python3

import csv
import math
import os
import pandas as pa
import re
import shutil
import six

import ancestor_go_extraction
import pathway_extraction.uniprot_retrieval_data as uniprot_retrieval_data

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'

class FileManagement():

    def __init__(self, name_of_the_file, already_analyzed_t_f):
        self._file_name, self._file_extension = os.path.splitext(name_of_the_file)
        self._already_analyzed_file_tf = self.string_to_boolean(already_analyzed_t_f)

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

    @property
    def already_analyzed_file_tf(self):
        return self._already_analyzed_file_tf

    @already_analyzed_file_tf.setter
    def already_analyzed_file_tf(self, boolean_response):
        self._already_analyzed_file_tf = self.string_to_boolean(boolean_response)

    def string_to_boolean(self, value):
        return value.lower() in ("yes", "true", "y", "t", "1")

    def go_label_number_dictionnary_creation(self, file_name, specification):
        d_go_label_to_number = {}

        with open (file_name, 'r') as file:
            query_results_file = file.read()
            query_results_file = query_results_file.replace(" ,\n", "\n")
            query_results_file = query_results_file.replace(" , ", "\t")
            query_results_modified = open(temporary_directory + "queryResults.tsv", "w")
            query_results_modified.write(query_results_file)
            query_results_modified.close()

        query_results_dataframe = pa.read_csv(temporary_directory + "queryResults.tsv", sep = "\t")

        query_results_dataframe.columns = [["subject", "label", "NarrowSynonym", "BroadSynonym", "RelatedSynonym"]]

        quote_deletion = lambda x: x.replace('"', '')
        query_results_dataframe["subject"] = query_results_dataframe["subject"].apply(quote_deletion)

        go_isolation = lambda x: x[32:]
        query_results_dataframe["subject"] = query_results_dataframe["subject"].apply(go_isolation)

        query_results_dataframe["subject"] = query_results_dataframe["subject"].str.replace("_", ":")

        if specification == "inverse":
            d_go_label_to_number = dict(zip(query_results_dataframe['subject'], query_results_dataframe['label']))

            return d_go_label_to_number

        else:
            space_deletion = lambda x: x.replace(" ", "")
            query_results_dataframe["label"] = query_results_dataframe["label"].apply(space_deletion)

            space_deletion = lambda x: str(x).replace(" ", "")
            query_results_dataframe["NarrowSynonym"] = query_results_dataframe["NarrowSynonym"].apply(space_deletion)
            query_results_dataframe["BroadSynonym"] = query_results_dataframe["BroadSynonym"].apply(space_deletion)
            query_results_dataframe["RelatedSynonym"] = query_results_dataframe["RelatedSynonym"].apply(space_deletion)

            d_go_label_to_number = dict(zip(query_results_dataframe['label'], query_results_dataframe['subject']))

        d_go_label_with_synonym = {}

        query_results_dataframe.set_index(query_results_dataframe['subject'], inplace = True)

        for index, row in query_results_dataframe.iterrows():
            if row['NarrowSynonym'] not in d_go_label_with_synonym and row['BroadSynonym'] not in d_go_label_with_synonym\
            and row['RelatedSynonym'] not in d_go_label_with_synonym:
                if row['NarrowSynonym'] != 'nan':
                    d_go_label_with_synonym[row['NarrowSynonym']] = index
                if row['BroadSynonym'] != 'nan':
                    d_go_label_with_synonym[row['BroadSynonym']] = index
                if row['RelatedSynonym'] != 'nan':
                    d_go_label_with_synonym[row['RelatedSynonym']] = index

        keys_to_delete = []

        for key, values in d_go_label_with_synonym.items():
            if values == []:
                keys_to_delete.append(key)

        for key in keys_to_delete:
            del d_go_label_with_synonym[key]

        return d_go_label_to_number, d_go_label_with_synonym

    def translate_go_label_into_go_number(self, gos_labels, d_go_label_to_number):
        gos_numbers = []

        if type(gos_labels) is float:
            if math.isnan(gos_labels):
                gos_labels = [str(gos_labels)]

        gos_labels = [go.replace(" ", "") for go in gos_labels]

        for go_label in gos_labels:
            if go_label in d_go_label_to_number:
                go_number = d_go_label_to_number[go_label]
                gos_numbers.append(go_number)
            else:
                gos_numbers.append(go_label)

        return gos_numbers

    def fix_obsolete_go_term(self, gos_labels_and_numbers, d_go_label_to_number, d_go_label_with_synonym):
        gos_numbers = []

        gos_labels_and_numbers = [go.replace(" ", "") for go in gos_labels_and_numbers]

        for go_label_and_number in gos_labels_and_numbers:
            if go_label_and_number == "-":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] == "GO":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] != 'GO' and go_label_and_number[0:2] != "-":
                go_labelObsolete = "obsolete" + go_label_and_number
                if go_labelObsolete in d_go_label_to_number:
                    go_number = d_go_label_to_number[go_labelObsolete]
                    gos_numbers.append(go_number)
                elif go_labelObsolete in d_go_label_with_synonym:
                    go_number = d_go_label_with_synonym[go_label_fixed]
                    gos_numbers.append(go_number)
                else:
                    gos_numbers.append(go_label_and_number)

        return gos_numbers

    def fix_wrong_n_or_l_term(self, gos_labels_and_numbers, d_go_label_to_number, d_go_label_with_synonym):
        gos_numbers = []

        gos_labels_and_numbers = [go.replace(" ", "") for go in gos_labels_and_numbers]

        for go_label_and_number in gos_labels_and_numbers:
            if go_label_and_number == "-":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] == "GO":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] != 'GO' and go_label_and_number[0:2] != "-":
                if "N" in go_label_and_number :
                    go_label_n_wrong = go_label_and_number.replace("N", "L")
                    if go_label_n_wrong in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_n_wrong]
                        gos_numbers.append(go_number)
                    elif go_label_n_wrong in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                    else:
                        gos_numbers.append(go_label_and_number)

                elif "L" in go_label_and_number :
                    go_label_l_wrong = go_label_and_number.replace("L", "N")
                    if go_label_l_wrong in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_l_wrong]
                        gos_numbers.append(go_number)
                    elif go_label_l_wrong in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)

                    else:
                        gos_numbers.append(go_label_and_number)
                else:
                    gos_numbers.append(go_label_and_number)

        return gos_numbers

    def fix_problems_with_synonym(self, gos_labels_and_numbers, d_go_label_to_number, d_go_label_with_synonym):
        gos_numbers = []

        gos_labels_and_numbers = [go.replace(" ", "") for go in gos_labels_and_numbers]

        for go_label_and_number in gos_labels_and_numbers:
            if go_label_and_number == "-":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] == "GO":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] != 'GO' and go_label_and_number[0:2] != "-":
                if go_label_and_number in d_go_label_with_synonym :
                        go_number = d_go_label_with_synonym[go_label_and_number]
                        gos_numbers.append(go_number)
                else:
                    gos_numbers.append(go_label_and_number)

        return gos_numbers

    def fix_terms_issue(self, gos_labels_and_numbers, d_go_label_to_number, d_go_label_with_synonym):
        gos_numbers = []

        gos_labels_and_numbers = [go.replace(" ", "") for go in gos_labels_and_numbers]

        for go_label_and_number in gos_labels_and_numbers:
            if go_label_and_number == "-":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] == "GO":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] != 'GO' and go_label_and_number[0:2] != "-":
                if "hydrogen" in go_label_and_number :
                    go_label_fixed = go_label_and_number.replace("hydrogen", "proton")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "al" in go_label_and_number :
                    go_label_fixed = go_label_and_number.replace("al", "")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "dependent" in go_label_and_number :
                    go_label_fixed = go_label_and_number.replace("dependent", "templated")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "sequence-specific" in go_label_and_number :
                    go_label_fixed = go_label_and_number.replace("sequence-specific", "") + ",sequence-specific"
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "homocysteineS-methyltransferaseactivity" in go_label_and_number:
                    go_label_fixed = "S-adenosylmethionine-" + go_label_and_number
                    if go_label_fixed in d_go_label_to_number[go_label_fixed]:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "Cul4-RINGubiquitinligasecomplex" in go_label_and_number:
                    go_label_fixed = go_label_and_number[:9] + "E3" + go_label_and_number[9:]
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "organ" in go_label_and_number:
                    go_label_fixed = "animal" + go_label_and_number
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "mitosis" in go_label_and_number:
                    go_label_fixed = "mitoticnucleardivision"
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "stimulus" in go_label_and_number:
                    go_label_fixed = go_label_and_number[:-len("stimulus")]
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "tailtipmorphogenesis" in go_label_and_number:
                    go_label_fixed = "nematodemale" + go_label_and_number
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "carboxylesteraseactivity" in go_label_and_number:
                    go_label_fixed = go_label_and_number[0:len('carboxyl')] + 'icesterhydrolaseactivity'
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "homophiliccelladhesion" in go_label_and_number:
                    go_label_fixed = go_label_and_number + 'viaplasmamembraneadhesionmolecules'
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "LSUrRNAbinding" in go_label_and_number:
                    go_label_fixed = 'largeribosomalsubunit' + go_label_and_number[0:-len('LSU')]
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "threonylcarbamoyladenosine" in go_label_and_number:
                    go_label_fixed = 'cyclic' + go_label_and_number
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "intraflagellartransportparticleB" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("flagellar", "ciliary")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "intraflagellartransportparticleA" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("flagellar", "ciliary")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "extracellularvesicularexosome" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("vesicular", "")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "Mphaseofmitoticcellcycle" in go_label_and_number:
                    go_label_fixed = 'mitoticMphase'
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "ATADPantiporteractivity" in go_label_and_number:
                    go_label_fixed = go_label_and_number[0:len('AT')] + 'P:' + go_label_and_number[len('AT'):]
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "ribonucleaseHactivity" in go_label_and_number:
                    go_label_fixed = "RNA-DNAhybridribonucleaseactivity"
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "methylatedhistoneresiduebinding" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("residue", "")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "ciliumaxonemeassembly" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("axoneme", "")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "telomerictemplateRNAreversetranscriptaseactivity" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("ictemplate", "ase")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "ciliumaxoneme" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("cilium", "")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "autophagicvacuolemembrane" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("icvacuole", "osome")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "autophagicvacuoleassembly" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("icvacuole", "osome")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "methylenetetrahydrofolatereductase(NADPH)activity" in go_label_and_number:
                    go_label_fixed = go_label_and_number.replace("P", "(P)")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                elif "cytoplasmictransport" in go_label_and_number:
                    go_label_fixed = go_label_and_number + ",nursecelltooocyte"
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                else:
                    gos_numbers.append(go_label_and_number)

        return gos_numbers

    def fix_dash_in_excess(self, gos_labels_and_numbers, d_go_label_to_number, d_go_label_with_synonym):
        gos_numbers = []

        gos_labels_and_numbers = [go.replace(" ", "") for go in gos_labels_and_numbers]

        for go_label_and_number in gos_labels_and_numbers:
            if go_label_and_number == "-":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] == "GO":
                gos_numbers.append(go_label_and_number)
            elif go_label_and_number[0:2] != 'GO' and go_label_and_number[0:2] != "-":
                if "-" in go_label_and_number :
                    go_label_fixed = go_label_and_number.replace("-", "")
                    if go_label_fixed in d_go_label_to_number:
                        go_number = d_go_label_to_number[go_label_fixed]
                        gos_numbers.append(go_number)
                    elif go_label_fixed in d_go_label_with_synonym:
                        go_number = d_go_label_with_synonym[go_label_fixed]
                        gos_numbers.append(go_number)
                else:
                    gos_numbers.append(go_label_and_number)

        return gos_numbers

    def go_translation(results_dataframe):
        d_go_label_to_number, d_go_label_with_synonym = self.go_label_number_dictionnary_creation(input_directory + "queryResults.csv", "normal")

        translation = lambda x: self.translate_go_label_into_go_number(x, d_go_label_to_number)
        results_dataframe['GOs'] = results_dataframe['GOs'].apply(translation)

        correction_Synonym_Issues = lambda x : self.fix_problems_with_synonym(x,  d_go_label_to_number, d_go_label_with_synonym)
        results_dataframe['GOs'] = results_dataframe['GOs'].apply(correction_Synonym_Issues)

        correction_obsolete_go = lambda x: self.fix_obsolete_go_term(x, d_go_label_to_number, d_go_label_with_synonym)
        results_dataframe['GOs'] = results_dataframe['GOs'].apply(correction_obsolete_go)

        correction_n_or_l = lambda x: self.fix_wrong_n_or_l_term(x, d_go_label_to_number, d_go_label_with_synonym)
        results_dataframe['GOs'] = results_dataframe['GOs'].apply(correction_n_or_l)

        correction_terms_issue = lambda x : self.fix_terms_issue(x, d_go_label_to_number, d_go_label_with_synonym)
        results_dataframe['GOs'] = results_dataframe['GOs'].apply(correction_terms_issue)

        correction_dash_issue = lambda x : self.fix_dash_in_excess(x, d_go_label_to_number, d_go_label_with_synonym)
        results_dataframe['GOs'] = results_dataframe['GOs'].apply(correction_dash_issue)

        return results_dataframe

    def cleaning_value(self, dataframe, value):
        value_dataframe = dataframe[dataframe.GOs.str.match(value) == True]

        dataframe.set_index("Gene_Name", inplace = True)
        for index in value_dataframe['Gene_Name'].tolist():
            dataframe = dataframe.drop(index)
        dataframe.reset_index(inplace = True)

        return dataframe

    def cleaning_nan_value(self, dataframe, column):
        for index, row in dataframe.iterrows():
            if type(row[column]) is float:
                if math.isnan(dataframe.get_value(index, column)):
                    dataframe = dataframe.drop([index])

        return dataframe

    def rewriting_file(self, newtable, file_name):
        newtable.to_csv(temporary_directory + file_name, "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

    def find_column_of_interest(self, df):
        columns = df.columns.tolist()

        go_label_expression = r"[FPC]{1}:[\w]*"
        go_number_expression = r"[FPC]{1}:GO[:_][\d]{7}"
        ec_expression = r"EC:[\d]{1}[\.]{1}[\d]{,2}[\.]{,1}[\d]{,2}[\.]{,1}[\d]{,2}"
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
            go_number_column = max(go_number_columns, key = go_number_columns.get)
            go_column = go_number_column
            go_label_columns.pop(go_number_column, None)
        ec_column = max(ec_columns, key = ec_columns.get)
        ipr_column = max(ipr_columns, key = ipr_columns.get)
        if ko_keggs:
            ko_kegg = max(ko_keggs, key = ko_keggs.get)

        go_label_column = max(go_label_columns, key = go_label_columns.get)

        if not go_number_columns:
            go_column = go_label_column

        return go_column, ec_column, ipr_column

    def column_go_cleaning(self, type_file):

        name_input_file = self.file_name
        extension_input_file = self.file_extension

        if extension_input_file == '.xls':
            results_dataframe = pa.read_excel(input_directory + name_input_file + extension_input_file, sep = None, na_values = "")
        else:
            results_dataframe = pa.read_csv(input_directory + name_input_file + extension_input_file, sep = None, engine = "python", na_values = "")

        yes_or_no = input("Is the first columns of your file, the column containing gene name? ")

        yes_answers = ['yes', 'y', 'oui', 'o']
        if yes_or_no.lower() in yes_answers :
            name_gene_column = results_dataframe.columns[0]
        else :
            name_gene_column = input("Write the name of the column containing the gene names : ")

        if type_file == "gene_list":
            results_dataframe = results_dataframe[[name_gene_column]]
            results_dataframe.columns = [['Gene_Name']]
            self.rewriting_file(results_dataframe, name_input_file + "GOsTranslatedAndFixed.tsv")
            return name_input_file + "GOsTranslatedAndFixed.tsv"

        go_column, ec_column, ipr_column = self.find_column_of_interest(results_dataframe)

        yes_or_no = input("Do you want to try to retrieve data from blast results? ").lower()
        if yes_or_no in yes_answers :
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

        if yes_or_no in yes_answers :
            results_dataframe = uniprot_retrieval_data.extract_information_from_uniprot(results_dataframe)

        results_dataframe = self.cleaning_value(results_dataframe, '-')
        results_dataframe = self.cleaning_nan_value(results_dataframe, 'GOs')

        self.rewriting_file(results_dataframe, name_input_file + "GOsTranslatedAndFixed.tsv")

        return name_input_file + "GOsTranslatedAndFixed.tsv"

class FileManagementGeneGOs(FileManagement):

    def __init__(self, name_of_the_file, already_analyzed_tf, type_of_the_file, column_name):
        FileManagement.__init__(self, name_of_the_file, already_analyzed_tf)
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

        df = pa.read_csv(temporary_directory + file_name_temporary, '\t')

        for index, row in df.iterrows():
            df[column_analyzed_object].loc[index] = ancestor_go_extraction.union_go_and_their_ancestor(row[column_analyzed_object].split(","))

        df.to_csv(temporary_directory + file_name_temporary, '\t', index = False)

    def create_gene_object_analysis_file(self, file_name, columns_names, column_analyzed_object):
        df = pa.read_csv(temporary_directory + file_name, sep = "\t")
        df = df[columns_names]
        df.set_index(columns_names[0], inplace = True)

        csvfile = open(temporary_directory + "Gene_" + file_name + column_analyzed_object + ".tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow((columns_names[0], columns_names[1]))

        for gene, list_of_list_of_analyzed_objects in df.iterrows():
            for list_of_analyzed_objects in list_of_list_of_analyzed_objects:
                for analyzed_object in list_of_analyzed_objects.split(","):
                    writer.writerow((gene, analyzed_object))

        csvfile.close()

    def file_gene_gos_gestion(self):
        analyzed_object_name = self.analyzed_object_name
        name_of_the_file = self.file_name
        extension_of_the_file  = self.file_extension
        type_of_the_file = self.type_file

        if self.already_analyzed_file_tf == True and type_of_the_file == 'genome':
            shutil.copy(input_directory + name_of_the_file + extension_of_the_file, temporary_directory) 
            file_name_temporary = name_of_the_file + extension_of_the_file
        if self.already_analyzed_file_tf == False:
            file_name_temporary = self.column_go_cleaning(type_of_the_file)

        if type_of_the_file == 'genome':
            if self.already_analyzed_file_tf == False:
                self.go_ancestors_list_of_interest(analyzed_object_name, file_name_temporary)
            counting_object_file = self.counting_genome(file_name_temporary, 'CountsReference', analyzed_object_name)

            return file_name_temporary, counting_object_file

        if self.type_file == 'gene_list':
            counting_object_file, number_of_gene = self.counting_gene_list(file_name_temporary, 'Counts', analyzed_object_name)

            return counting_object_file, number_of_gene

class FileManagementGeneGOsGenome(FileManagementGeneGOs):

    def __init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name):
        FileManagementGeneGOs.__init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name)

    def counting_genome(self, file_name_temporary, column_name, column_analyzed_object):
        analyzed_objects = []
        df = pa.read_csv(temporary_directory + file_name_temporary, sep="\t")

        for index, row in df.iterrows():
            for analyzed_object in row[column_analyzed_object].split(","):
                analyzed_objects.append(analyzed_object)

        counts_df_Genome = pa.DataFrame(analyzed_objects)
        counts_df_Genome.columns = [column_analyzed_object]
        counts_df_Genome = counts_df_Genome.groupby(column_analyzed_object).size().rename(column_name)
        counts_df_Genome = counts_df_Genome.to_frame()

        counts_df_Genome.reset_index(inplace = True)
        self.rewriting_file(counts_df_Genome, "counting_objects_in_genome.tsv")

        return "counting_objects_in_genome"

class FileManagementGeneGOsInterest(FileManagementGeneGOs):

    def __init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name, file_name_genome):
        FileManagementGeneGOs.__init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name)
        self._file_genome_reference_name = file_name_genome

    @property
    def genome_file_reference_name(self):
        return self._file_genome_reference_name

    @genome_file_reference_name.setter
    def genome_file_reference_name(self, file_name):
        self._file_genome_reference_name = file_name

    def counting_gene_list(self, file_name_temporary, column_name, column_analyzed_object):
        analyzed_objects = []
        df = pa.read_csv(temporary_directory + file_name_temporary, sep = "\t")
        df = df[['Gene_Name']]
        df.set_index("Gene_Name", inplace = True)

        df_genome = pa.read_csv(temporary_directory + self.genome_file_reference_name, sep = "\t")
        df_genome = df_genome[['Gene_Name', 'GOs']]
        df_genome.set_index("Gene_Name", inplace = True)

        df_joined = df.join(df_genome)
        df_joined = df_joined[pa.notnull(df_joined[column_analyzed_object])]
        df_joined.reset_index(inplace = True)
        self.rewriting_file(df_joined, file_name_temporary)
        df_joined.set_index("Gene_Name", inplace = True)

        for index, row in df_joined.iterrows():
            for analyzed_object in row[column_analyzed_object].split(","):
                analyzed_objects.append(analyzed_object)

        counts_df = pa.DataFrame(analyzed_objects)
        counts_df.columns = [column_analyzed_object]
        counts_df = counts_df.groupby(column_analyzed_object).size().rename(column_name)
        counts_df = counts_df.to_frame()

        numberOfGene = len(df.index.unique())

        counts_df.reset_index(inplace = True)
        self.rewriting_file(counts_df, "counting_objects_in_interest.tsv")

        return "counting_objects_in_interest", numberOfGene

class FileManagementGeneGO(FileManagement):

    def __init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name):
        FileManagement.__init__(self, name_of_the_file, already_analyzed_t_f)
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

    def file_gene_gos_gestion():
        analyzed_object_name = self.analyzed_object_name
        name_of_file = self.file_name

        self.column_go_cleaning()

        if self.type_file == 'gene_list':
            self.genome_file_processing(name_of_file)
        if self.type_file == 'genome':
            self.file_of_interest_processing(name_of_file)

class FileManagementGeneGOGenome(FileManagementGeneGO):

    def __init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name):
        FileManagementGeneGO.__init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name)

    def genome_file_processing(self, genome_file_name):
        df = pa.read_csv(temporary_directory + self.file_name + self.file_extension, sep = "\t", header = None)
        df.columns = [['Gene', 'GOs']]

        df .set_index("Gene", inplace = True)

        genes_gos_ancestors = {}
        for gene, row in df.iterrows():
            go_with_ancestor = ancestor_go_extraction.union_go_and_their_ancestor([row['GOs']])
            if gene not in genes_gos_ancestors:
                genes_gos_ancestors[gene] = go_with_ancestor
            else:
                genes_gos_ancestors[gene].extend(go_with_ancestor)

        csvfile = open(temporary_directory + self.file_name + "_with_ancestor.tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["Gene", "GO"])

        for gene in genes_gos_ancestors:
            gene_go_unique = []
            for go in genes_gos_ancestors[gene]:
                if go not in gene_go_unique:
                    writer.writerow([gene, go])
                    gene_go_unique.append(go)

        csvfile.close()

        df = pa.read_csv(temporary_directory + self.file_name + "_with_ancestor.tsv", sep = "\t", header = None)
        df.columns = [['Gene', 'GOs']]
        df.set_index("Gene", inplace = True)
        go_counts = {}

        for gene, row in df.iterrows():
            if row["GOs"] not in go_counts:
                go_counts[row["GOs"]] = 1
            elif row["GOs"] in go_counts:
                go_counts[row["GOs"]] += 1

        csvfile = open(temporary_directory + "number_go_in_genome.tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["GOs", "CountsReference"])

        for go in go_counts:
            writer.writerow([go, go_counts[go]])

        csvfile.close()

class FileManagementGeneGOInterest(FileManagementGeneGO):

    def __init__(self, name_of_the_file, already_analyzed_t_f, type_of_the_file, column_name, file_name_genome):
        FileManagement.__init__(self, name_of_the_file, already_analyzed_t_f)
        self._file_genome_reference_name = file_name_genome

    @property
    def genome_file_reference_name(self):
        return self._file_genome_reference_name

    @genome_file_reference_name.setter
    def genome_file_reference_name(self, file_name):
        self._file_genome_reference_name = file_name

    def file_of_interest_processing(self):
        df = pa.read_csv(temporary_directory + self.file_name + self.file_extension, sep = "\t", header = None)
        df.columns = [['Gene']]
        df.set_index("Gene", inplace = True)

        df_genome = pa.read_csv(temporary_directory + self.genome_file_reference_name + "_with_ancestor.tsv", sep = "\t", header = None)
        df_genome.columns = [['Gene', 'GOs']]
        df_genome.set_index("Gene", inplace = True)

        df_joined = df.join(df_genome)
        go_counts = {}

        for gene, row in df_joined.iterrows():
            if row["GOs"] not in go_counts:
                go_counts[row["GOs"]] = 1
            elif row["GOs"] in go_counts:
                go_counts[row["GOs"]] += 1

        csvfile = open(temporary_directory + "number_go_in_interest.tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["GOs", "Counts"])

        for go in go_counts:
            writer.writerow([go, go_counts[go]])

        csvfile.close()

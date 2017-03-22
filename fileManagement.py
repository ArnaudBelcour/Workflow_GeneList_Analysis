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

    def dict_to_file(dictionary, file_name):
        df = pa.DataFrame.from_dict(dictionary, orient='index')
        df.reset_index(inplace=True)
        df.columns = [['GOlabel', 'GOnumber']]
        df.to_csv(temporary_directory + file_name + '.tsv', sep='\t', index=False)

    def go_label_number_dictionnary_creation_from_http(self, specification='normal'):
        '''
            Create a dictionnary containing GO labels (as key) associated with their GO numbers (as value), if specification is 'normal'. Use to translate GO labels into GO numbers.
            Or GO numbers (as key) associated with their GO labels (as value), if specification is 'inverse'. Use to translate GO numbers into GO labels.
            Create also a dictionnary containing the synonym of the GO terms.
        '''
        d_go_label_to_number = {}

        go_ontology = pronto.Ontology('http://purl.obolibrary.org/obo/go/go-basic.obo')

        if specification == "inverse":
            for go_term in go_ontology:
                d_go_label_to_number[go_term.id] = go_term.name
            dict_to_file(d_go_label_to_number, 'go_number_label')

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

    def go_label_number_dictionnary_creation_from_file(self, specification='normal'):
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
        d_go_label_to_number, d_go_label_with_synonym = self.go_label_number_dictionnary_creation()

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
        ec_column = max(ec_columns, key=ec_columns.get)
        ipr_column = max(ipr_columns, key=ipr_columns.get)
        if ko_keggs:
            ko_kegg = max(ko_keggs, key=ko_keggs.get)

        if not go_number_columns:
            go_label_column = max(go_label_columns, key=go_label_columns.get)
            go_column = go_label_column

        return go_column, ec_column, ipr_column

    def column_data_cleaning(self, type_file):

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
            results_dataframe = pa.read_csv(input_directory + name_input_file + extension_input_file, sep=None, engine="python", na_values="")

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

        results_dataframe.to_csv(temporary_directory + name_input_file + "GOsTranslatedAndFixed.tsv", "\t", index=True, quoting=csv.QUOTE_NONE)

        return name_input_file + "GOsTranslatedAndFixed.tsv"

class FileManagementGeneGOs(FileManagement):

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

    def create_gene_object_analysis_file(self, file_name, columns_names, column_analyzed_object):
        df = pa.read_csv(temporary_directory + file_name, sep="\t")
        df = df[columns_names]
        df.set_index(columns_names[0], inplace=True)

        csvfile = open(temporary_directory + "Gene_" + file_name + column_analyzed_object + ".tsv", "w", newline="")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow((columns_names[0], columns_names[1]))

        for gene, list_of_list_of_analyzed_objects in df.iterrows():
            for list_of_analyzed_objects in list_of_list_of_analyzed_objects:
                for analyzed_object in list_of_analyzed_objects.split(","):
                    writer.writerow((gene, analyzed_object))

        csvfile.close()

    def file_gene_gos_gestion(self):
        string_to_boolean = lambda value: value.lower() in ("yes", "true", "y", "t", "1")

        analyzed_object_name = self.analyzed_object_name
        name_of_the_file = self.file_name
        extension_of_the_file  = self.file_extension
        type_of_the_file = self.type_file

        if type_of_the_file == 'genome':
            already_analyzed_file_yes_no = string_to_boolean(input("Does this file have already been analyzed? "))
            if already_analyzed_file_yes_no == True:
                shutil.copy(input_directory + name_of_the_file + extension_of_the_file, temporary_directory)
                file_name_temporary = name_of_the_file + extension_of_the_file

            elif already_analyzed_file_yes_no == False:
                file_name_temporary = self.column_data_cleaning(type_of_the_file)
                self.go_ancestors_list_of_interest(analyzed_object_name, file_name_temporary)

                session = requests.Session()
                pathway_extractor.data_retrieval_from_GO(file_name_temporary)
                pathway_extractor.main(file_name_temporary, session)

                mapping_pathway_data.main(file_name_temporary)

            counting_object_file, number_of_gene_genome = self.counting_genome(file_name_temporary, 'CountsReference', analyzed_object_name)

            return file_name_temporary, counting_object_file, number_of_gene_genome

        elif type_of_the_file == 'gene_list':
            file_name_temporary = self.column_data_cleaning(type_of_the_file)
            counting_object_file, number_of_gene_list = self.counting_gene_list(file_name_temporary, 'Counts', analyzed_object_name)

            return counting_object_file, number_of_gene_list

class FileManagementGeneGOsGenome(FileManagementGeneGOs):

    def __init__(self, name_of_the_file, type_of_the_file, column_name):
        FileManagementGeneGOs.__init__(self, name_of_the_file, type_of_the_file, column_name)

    def counting_genome(self, file_name_temporary, column_name, column_analyzed_object):
        analyzed_objects = []
        df = pa.read_csv(temporary_directory + file_name_temporary, sep="\t")
        df.replace(np.nan, '', regex=True, inplace=True)
        df.set_index('Gene_Name', inplace=True)

        #for index, row in df.iterrows():
            #if row['GOs'] == '':
                #df.drop(index, inplace=True)

        for index, row in df.iterrows():
            for analyzed_object in row[column_analyzed_object].split(","):
                analyzed_objects.append(analyzed_object)

        counts_df_genome = pa.DataFrame(analyzed_objects)
        counts_df_genome.columns = [column_analyzed_object]
        counts_df_genome = counts_df_genome.groupby(column_analyzed_object).size().rename(column_name)
        counts_df_genome = counts_df_genome.to_frame()
        #counts_df_genome['CountsProportion'] = counts_df_genome['CountsReference']/len(df.index)
        counts_df_genome.reset_index(inplace=True)
        counts_df_genome.to_csv(temporary_directory + "counting_objects_in_genome.tsv", "\t", index=False, header=True, quoting=csv.QUOTE_NONE)

        number_of_gene = len(df.index.unique())

        return "counting_objects_in_genome", number_of_gene

class FileManagementGeneGOsInterest(FileManagementGeneGOs):

    def __init__(self, name_of_the_file, type_of_the_file, column_name, file_name_genome):
        FileManagementGeneGOs.__init__(self, name_of_the_file, type_of_the_file, column_name)
        self._file_genome_reference_name = file_name_genome

    @property
    def genome_file_reference_name(self):
        return self._file_genome_reference_name

    @genome_file_reference_name.setter
    def genome_file_reference_name(self, file_name):
        self._file_genome_reference_name = file_name

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

        #for index, row in df_joined.iterrows():
            #if row['GOs'] == '':
                #df_joined.drop(index, inplace=True)

        for index, row in df_joined.iterrows():
            for analyzed_object in row[column_analyzed_object].split(","):
                analyzed_objects.append(analyzed_object)

        counts_df = pa.DataFrame(analyzed_objects)
        counts_df.columns = [column_analyzed_object]
        counts_df = counts_df.groupby(column_analyzed_object).size().rename(column_name)
        counts_df = counts_df.to_frame()
        #counts_df['CountsProportion'] = counts_df['Counts']/len(df_joined.index)
        #counts_ref = pa.read_csv(temporary_directory + "counting_objects_in_genome.tsv", sep='\t')
        #counts_ref.set_index('GOs', inplace=True)
        #counts_df['Counts'] = counts_df['Counts'] - (56 * counts_df['CountsProportion'])
        number_of_gene = len(df.index.unique())

        counts_df.reset_index(inplace=True)
        counts_df.to_csv(temporary_directory + "counting_objects_in_interest.tsv", "\t", index=False, header=True, quoting=csv.QUOTE_NONE)

        return "counting_objects_in_interest", number_of_gene

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

    def file_gene_gos_gestion():
        analyzed_object_name = self.analyzed_object_name
        name_of_file = self.file_name

        self.column_data_cleaning()

        if self.type_file == 'gene_list':
            self.genome_file_processing(name_of_file)
        if self.type_file == 'genome':
            self.file_of_interest_processing(name_of_file)

class FileManagementGeneGOGenome(FileManagementGeneGO):

    def __init__(self, name_of_the_file, type_of_the_file, column_name):
        FileManagementGeneGO.__init__(self, name_of_the_file, type_of_the_file, column_name)

    def genome_file_processing(self, genome_file_name):
        df = pa.read_csv(temporary_directory + self.file_name + self.file_extension, sep="\t", header=None)
        df.columns = [['Gene', 'GOs']]

        df .set_index("Gene", inplace=True)

        genes_gos_ancestors = {}
        for gene, row in df.iterrows():
            go_with_ancestor = ancestor_go_extraction.union_go_and_their_ancestor([row['GOs']], 'list')
            if gene not in genes_gos_ancestors:
                genes_gos_ancestors[gene] = go_with_ancestor
            else:
                genes_gos_ancestors[gene].extend(go_with_ancestor)

        csvfile = open(temporary_directory + self.file_name + "_with_ancestor.tsv", "w", newline="")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["Gene", "GO"])

        for gene in genes_gos_ancestors:
            gene_go_unique = []
            for go in genes_gos_ancestors[gene]:
                if go not in gene_go_unique:
                    writer.writerow([gene, go])
                    gene_go_unique.append(go)

        csvfile.close()

        df = pa.read_csv(temporary_directory + self.file_name + "_with_ancestor.tsv", sep = "\t", header=None)
        df.columns = [['Gene', 'GOs']]
        df.set_index("Gene", inplace=True)
        go_counts = {}

        for gene, row in df.iterrows():
            if row["GOs"] not in go_counts:
                go_counts[row["GOs"]] = 1
            elif row["GOs"] in go_counts:
                go_counts[row["GOs"]] += 1

        csvfile = open(temporary_directory + "counting_objects_in_genome.tsv", "w", newline="")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["GOs", "CountsReference"])

        for go in go_counts:
            writer.writerow([go, go_counts[go]])

        csvfile.close()

        return "counting_objects_in_genome"

class FileManagementGeneGOInterest(FileManagementGeneGO):

    def __init__(self, name_of_the_file, type_of_the_file, column_name, file_name_genome):
        FileManagement.__init__(self, name_of_the_file)
        self._file_genome_reference_name = file_name_genome

    @property
    def genome_file_reference_name(self):
        return self._file_genome_reference_name

    @genome_file_reference_name.setter
    def genome_file_reference_name(self, file_name):
        self._file_genome_reference_name = file_name

    def file_of_interest_processing(self):
        df = pa.read_csv(temporary_directory + self.file_name + self.file_extension, sep="\t", header=None)
        df.columns = [['Gene']]
        df.set_index("Gene", inplace=True)

        df_genome = pa.read_csv(temporary_directory + self.genome_file_reference_name + "_with_ancestor.tsv", sep="\t", header=None)
        df_genome.columns = [['Gene', 'GOs']]
        df_genome.set_index("Gene", inplace=True)

        df_joined = df.join(df_genome)
        go_counts = {}

        for gene, row in df_joined.iterrows():
            if row["GOs"] not in go_counts:
                go_counts[row["GOs"]] = 1
            elif row["GOs"] in go_counts:
                go_counts[row["GOs"]] += 1

        csvfile = open(temporary_directory + "counting_objects_in_interest.tsv", "w", newline="")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow(["GOs", "Counts"])

        for go in go_counts:
            writer.writerow([go, go_counts[go]])

        csvfile.close()

        return 'counting_objects_in_interest'

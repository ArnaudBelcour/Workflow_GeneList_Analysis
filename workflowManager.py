#!/usr/bin/env python3

import os
import six

from enrichmentAnalysis import GOEnrichmentAnalysis
from fileManagement import FileManagementGeneGOs

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
temporary_directory_database = '../temporaryFiles/databases/'
output_directory = 'outputFiles/'

def workflow_mainager():
    if os.path.exists(input_directory[:-1]) == False :
        os.makedirs(input_directory)
        sys.exit("No input data, please put your data fiels in inputFiles directory.")
    if os.path.exists(temporary_directory[:-1]) == False :
        os.makedirs(temporary_directory)
        os.makedirs(temporary_directory_database)
    if os.path.exists(output_directory[:-1]) == False :
        os.makedirs(output_directory)

    if not os.listdir(input_directory):
        sys.exit("No input data, please put your data fiels in inputFiles directory.")

    sentence_choice = "Do you want to download database datas? "
    yes_or_no = input(sentence_choice)
    yes_answers = ['yes', 'y', 'oui', 'o']

    if yes_or_no in yes_answers:
        input_file_of_interest_management.go_ancestors_list_of_interest(object_to_analyze)

    sentence_choice = "Write the name of your input file containing differentially expressed gene : "
    name_de_input_file = input(sentence_choice)

    sentence_choice = "Write the name of your input file containing genome : "
    name_reference_input_file = input(sentence_choice)

    input_file_gestion = FileManagementGeneGOs(name_de_input_file, name_reference_input_file, 'GOs')

    input_file_gestion.file_gene_gos_gestion()

    d_go_label_to_number = input_file_gestion.go_label_number_dictionnary_creation(input_directory + "queryResults.csv", 'inverse')

    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', d_go_label_to_number)
    object_to_analyze = go_enrichment_analysis.get_object_to_analyze()

    sentence_choice_number_gene = "Enter the number of genes in the genome of your organism : "
    number_of_genesInGenome = int(input(sentence_choice_number_gene))

    sentence_choice_alpha = "Enter the alpha risk : "
    alpha = float(input(sentence_choice_alpha))

    go_enrichment_analysis.set_file_of_interest(file_of_interest_name)
    go_enrichment_analysis.set_file_of_reference(fileOfGenomeName)
    go_enrichment_analysis.set_number_of_analyzed_object_of_interest(number_of_gene)
    go_enrichment_analysis.set_number_of_analyzed_object_of_reference(number_of_genesInGenome)
    go_enrichment_analysis.set_alpha(alpha)
    go_enrichment_analysis.enrichment_analysis()

workflow_mainager()

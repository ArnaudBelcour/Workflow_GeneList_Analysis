#!/usr/bin/env python3

import os
import six

from enrichmentAnalysis import GOEnrichmentAnalysis
from fileManagement import FileManagement

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
    nameDEInputFile = input(sentence_choice)
    input_file_of_interest_management = FileManagement(nameDEInputFile)
    file_of_interest_name_input = input_file_of_interest_management.get_file_name()

    d_go_label_to_number, d_go_label_with_synonym = input_file_of_interest_management.go_label_number_dictionnary_creation(input_directory + "queryResults.csv", 'inverse')

    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', d_go_label_to_number)
    object_to_analyze = go_enrichment_analysis.get_object_to_analyze()

    input_file_of_interest_management.column_go_cleaning()
    input_file_of_interest_management.create_gene_object_analysis_file(file_of_interest_name_input + object_to_analyze + "TranslatedAndFixed.tsv", ['Gene_Name', object_to_analyze], object_to_analyze)
    input_file_of_interest_management.go_ancestors_list_of_interest(object_to_analyze)

    sentence_choice = "Write the name of your input file containing genome : "
    name_reference_input_file = input(sentence_choice)
    input_reference_file_management = FileManagement(name_reference_input_file)
    file_referene_name_input = input_reference_file_management.get_file_name()

    input_reference_file_management.column_go_cleaning(go_enrichment_analysis)
    input_reference_file_management.create_gene_object_analysis_file(file_referene_name_input + object_to_analyze + "TranslatedAndFixed.tsv", ['Gene_Name', object_to_analyze], object_to_analyze)
    input_reference_file_management.go_ancestors_list_of_interest(object_to_analyze)

    file_of_interest_name, number_of_gene = input_file_of_interest_management.counting_gene_list(file_of_interest_name_input, 'Counts', object_to_analyze)
    fileOfGenomeName = input_reference_file_management.counting_genome(file_referene_name_input, 'CountsGenome', object_to_analyze)

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

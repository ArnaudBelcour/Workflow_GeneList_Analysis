#!/usr/bin/env python3

import os
import six

import go_extract_synonyms_label

from enrichmentAnalysis import GOEnrichmentAnalysis
from fileManagement import FileManagementGeneGOsGenome, FileManagementGeneGOsInterest

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'
output_directory = 'outputFiles/'

def workflow_mainager():
    if os.path.exists(input_directory[:-1]) == False :
        os.makedirs(input_directory)
        sys.exit("No input data, please put your data files in inputFiles directory.")
    if os.path.exists(temporary_directory[:-1]) == False :
        os.makedirs(temporary_directory)
        os.makedirs(temporary_directory_database)
    if os.path.exists(output_directory[:-1]) == False :
        os.makedirs(output_directory)

    if not os.listdir(input_directory):
        sys.exit("No input data, please put your data files in inputFiles directory.")

    yes_or_no = input("Do you want to recreate the GO terms and synonyms file? ")
    yes_answers = ['yes', 'y', 'oui', 'o']

    if yes_or_no in yes_answers:
        go_extract_synonyms_label.main()

    name_reference_input_file = input("Write the name of your input file containing genome : ")
    input_genome_file_gestion = FileManagementGeneGOsGenome(name_reference_input_file, 'genome' , 'GOs')
    reference_file_name, counting_reference_file_name, number_of_gene_genome = input_genome_file_gestion.file_gene_gos_gestion()

    name_de_input_file = input("Write the name of your input file containing differentially expressed gene : ")
    input_listde_file_gestion = FileManagementGeneGOsInterest(name_de_input_file, 'gene_list', 'GOs', reference_file_name)
    counting_interest_file_name, number_of_gene_list = input_listde_file_gestion.file_gene_gos_gestion()

    d_go_label_to_number = input_listde_file_gestion.go_label_number_dictionnary_creation('inverse')

    sentence_choice_alpha = "Enter the alpha risk : "
    alpha = float(input(sentence_choice_alpha))

    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', counting_interest_file_name, counting_reference_file_name, number_of_gene_list, number_of_gene_genome, alpha, 10000, d_go_label_to_number)
    go_enrichment_analysis.enrichment_analysis()

workflow_mainager()

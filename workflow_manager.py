#!/usr/bin/env python3

import os
import six

from enrichmentAnalysis import GOEnrichmentAnalysis
from fileManagement import FileManagementGeneGOGenome, FileManagementGeneInterest

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

    name_reference_input_file = input("Write the name of your input file containing genome : ")
    input_genome_file_gestion = FileManagementGeneGOGenome(name_reference_input_file, 'genome' , 'GOs')
    reference_file_name, counting_reference_file_name, number_of_gene_genome = input_genome_file_gestion.genome_file_processing()

    name_de_input_file = input("Write the name of your input file containing differentially expressed gene : ")
    input_listde_file_gestion = FileManagementGeneInterest(name_de_input_file, 'gene_list', 'GOs', reference_file_name)
    counting_interest_file_name, number_of_gene_list = input_listde_file_gestion.interest_file_processing()

    d_go_label_to_number = input_listde_file_gestion.go_label_number_dictionary_creation('inverse')

    sentence_choice_alpha = "Enter the alpha risk : "
    alpha = float(input(sentence_choice_alpha))

    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', counting_interest_file_name, counting_reference_file_name, number_of_gene_list, number_of_gene_genome, alpha, 10000, d_go_label_to_number)
    go_enrichment_analysis.enrichment_analysis()

workflow_mainager()

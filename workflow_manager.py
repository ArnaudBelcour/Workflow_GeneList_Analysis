#!/usr/bin/env python3

import os
import six

from enrichmentAnalysis import GOEnrichmentAnalysis
from fileManagement import FileManagementGeneGOsGenome, FileManagementGeneGOsInterest

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'
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

    yes_or_no = input("Do you want to download database datas? ")
    yes_answers = ['yes', 'y', 'oui', 'o']

    if yes_or_no in yes_answers:
        input_file_of_interest_management.go_ancestors_list_of_interest(object_to_analyze)

    name_reference_input_file = input("Write the name of your input file containing genome : ")
    already_analyzed_file_yes_no = input("Does this file already been analyzed? ")
    input_genome_file_gestion = FileManagementGeneGOsGenome(name_reference_input_file, already_analyzed_file_yes_no, 'genome' , 'GOs')
    reference_file_name, counting_reference_file_name = input_genome_file_gestion.file_gene_gos_gestion()

    name_de_input_file = input("Write the name of your input file containing differentially expressed gene : ")
    already_analyzed_file_yes_no = input("Does this file already been analyzed? ")
    input_listde_file_gestion = FileManagementGeneGOsInterest(name_de_input_file, already_analyzed_file_yes_no, 'gene_list', 'GOs', reference_file_name)
    counting_interest_file_name, number_of_gene = input_listde_file_gestion.file_gene_gos_gestion()

    d_go_label_to_number = input_listde_file_gestion.go_label_number_dictionnary_creation(input_directory + "queryResults.csv", 'inverse')

    sentence_choice_number_gene = "Enter the number of genes in the genome of your organism : "
    number_of_genes_in_genome = int(input(sentence_choice_number_gene))

    sentence_choice_alpha = "Enter the alpha risk : "
    alpha = float(input(sentence_choice_alpha))

    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', counting_interest_file_name, counting_reference_file_name, number_of_gene, number_of_genes_in_genome, alpha, 10000, d_go_label_to_number)
    go_enrichment_analysis.enrichment_analysis()

workflow_mainager()

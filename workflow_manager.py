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

    sentence_choice = "Do you want to download database datas? "
    yes_or_no = input(sentence_choice)
    yes_answers = ['yes', 'y', 'oui', 'o']

    if yes_or_no in yes_answers:
        input_file_of_interest_management.go_ancestors_list_of_interest(object_to_analyze)

    sentence_choice = "Write the name of your input file containing genome : "
    name_reference_input_file = input(sentence_choice)
    input_genome_file_gestion = FileManagementGeneGOsGenome(name_reference_input_file, 'genome' , 'GOs')
    file_of_genome_name = input_genome_file_gestion.file_gene_gos_gestion()

    sentence_choice = "Write the name of your input file containing differentially expressed gene : "
    name_de_input_file = input(sentence_choice)
    input_listde_file_gestion = FileManagementGeneGOsInterest(name_de_input_file, 'gene_list', 'GOs',
                                                                name_reference_input_file[:-len(input_genome_file_gestion.get_file_extension())])
    file_of_interest_name, number_of_gene = input_listde_file_gestion.file_gene_gos_gestion()

    d_go_label_to_number = input_listde_file_gestion.go_label_number_dictionnary_creation(input_directory + "queryResults.csv", 'inverse')

    sentence_choice_number_gene = "Enter the number of genes in the genome of your organism : "
    number_of_genes_in_genome = int(input(sentence_choice_number_gene))

    sentence_choice_alpha = "Enter the alpha risk : "
    alpha = float(input(sentence_choice_alpha))

    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', d_go_label_to_number, file_of_interest_name, file_of_genome_name, number_of_gene, number_of_genes_in_genome, alpha, 10000)
    go_enrichment_analysis.enrichment_analysis()

workflow_mainager()

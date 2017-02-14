#!/usr/bin/env python3

import csv
import pandas as pa
import re

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'

def ghost_koala_file_gestion(file_name):
    title_expression = r'\s+([A-Z]{1}[\D]+$)'
    pathway_expression = r'\s+[\d]{5}\s{1}'
    ko_expression = r'\s+K[\d]{5}'
    gene_expression = r'\s+Pldbra_eH_r[\d]{1}s[\d]{3}g[\d]{5}'

    title_pathways = {}
    pathway_kos = {}
    ko_genes = {}
    genes = []

    with open(input_directory + file_name) as file:
        for line in file:
            if re.match(gene_expression, line):
                genes_extracted = line.strip().split(",")
                for gene in genes_extracted:
                    genes.append(gene)
            if re.match(ko_expression, line):
                ko = line.strip()
                genes = []
                ko_genes[ko] = genes
            if re.match(pathway_expression, line):
                pathway = line.strip()
                ko_genes = {}
                pathway_kos[pathway] = ko_genes
            if re.match(title_expression, line):
                title = line.strip()
                pathway_kos = {}
                title_pathways[title] = pathway_kos

    title_pathways['Development']

    csvfile = open(temporary_directory + "ghost_koala.tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['Title', 'Pathway', 'KO', 'Gene'])

    for title in title_pathways:
        for pathway in title_pathways[title]:
            for ko in title_pathways[title][pathway]:
                for gene in title_pathways[title][pathway][ko]:
                    writer.writerow([title, pathway, ko, gene])

    csvfile.close()

def file_cleaning():
    df = pa.read_csv(temporary_directory + "ghost_koala.tsv", sep = "\t")

    series_KO = df.groupby('Gene')['KO'].apply(list)
    series_pathway = df.groupby('Gene')['Pathway'].apply(list)
    series_title = df.groupby('Gene')['Title'].apply(list)

    df_KO = series_KO.to_frame()
    df_pathway = series_pathway.to_frame()
    df_title = series_title.to_frame()

    df_KO['kegg_pathway'] = df_pathway['Pathway']
    df_KO['kegg_title'] = df_title['Title']

    df_KO['kegg_pathway'] = [['path:map' + pathway[:5] for pathway in pathways] for pathways in df_KO['kegg_pathway']]
    df_KO.to_csv(temporary_directory_database + "gene_with_kegg_pathway.tsv", sep = "\t")

ghost_koala_file_gestion('gene_pathwayq_GhostKoala.txt')
file_cleaning()
#!/usr/bin/env python3

import csv
import pandas as pa
import re

from tqdm import tqdm

from . import *

def ghost_koala_file_gestion(file_name):
    '''
        This script takes a txt file correspind to the extended version of a result of a Ghost Koala analysis.
        It uses the differents regular expressions to retrieves pathway associated with gene.
        To store the data a dictionnary of dictionnary of dictionnary of list is used.
        Genes names are stored in the list.
        K0XXXX associated with the gene are stored in the dictionnary as key and the value are the gene list.
        pathway ID are associated with the KO in the upper dictionnary as key and the value are K0XXXX.
        The title corresponds to the name of the generic process and is stored in a dictionnary as key and the value is the pathway ID.
    '''
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
            if re.match(title_expression, line):
                title = line.strip()
                pathway_kos = {}
                title_pathways[title] = pathway_kos
            if re.match(pathway_expression, line):
                pathway = line.strip()
                ko_genes = {}
                pathway_kos[pathway] = ko_genes
            if re.match(ko_expression, line):
                ko = line.strip()
                genes = []
                ko_genes[ko] = genes
            if re.match(gene_expression, line):
                genes_extracted = line.strip().split(",")
                for gene in genes_extracted:
                    genes.append(gene)

    csvfile = open(temporary_directory + "ghost_koala.tsv", "w", newline="")
    writer = csv.writer(csvfile, delimiter='\t')
    writer.writerow(['Title', 'Pathway', 'KO', 'Gene'])

    for title in tqdm(title_pathways):
        for pathway in title_pathways[title]:
            for ko in title_pathways[title][pathway]:
                for gene in title_pathways[title][pathway][ko]:
                    writer.writerow([title, pathway, ko, gene])

    csvfile.close()

def file_cleaning():
    df = pa.read_csv(temporary_directory + "ghost_koala.tsv", sep='\t')

    series_KO = df.groupby('Gene')['KO'].apply(list)
    series_pathway = df.groupby('Gene')['Pathway'].apply(list)
    series_title = df.groupby('Gene')['Title'].apply(list)

    df_KO = series_KO.to_frame()
    df_pathway = series_pathway.to_frame()
    df_title = series_title.to_frame()

    df_KO['kegg_pathway'] = df_pathway['Pathway']
    df_KO['kegg_title'] = df_title['Title']

    df_KO['kegg_pathway'] = [['path:map' + pathway[:5] for pathway in pathways] for pathways in df_KO['kegg_pathway']]
    df_KO.to_csv(temporary_directory_database + "gene_with_kegg_pathway.tsv", sep='\t')

def main():
    ghost_koala_file_gestion('gene_pathwayq_GhostKoala.txt')
    file_cleaning()

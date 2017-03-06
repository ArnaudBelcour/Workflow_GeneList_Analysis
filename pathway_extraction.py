#!/usr/bin/env python3

import os
import pandas as pa
import re
import subprocess

import pathway_extraction.chebi_from_go as chebi_from_go
import pathway_extraction.database_mapping_from_gos as database_mapping_from_gos
import pathway_extraction.eupathdb_pathway_extraction as eupathdb_pathway_extraction
import pathway_extraction.ghost_koala_pathway_extraction as ghost_koala_pathway_extraction
import pathway_extraction.interpro_pathway_extraction as interpro_pathway_extraction
import pathway_extraction.panther_pathway_mapping_uniprot as panther_pathway_mapping_uniprot
import pathway_extraction.reactome_pathway_extraction as reactome_pathway_extraction
import pathway_extraction.sparql_query_reactome_pathway_name as sparql_query_reactome_pathway_name

temporary_directory = 'temporaryFiles/'
temporary_directory_database = 'temporaryFiles/databases/'

def ipr_extraction(df_genome):
    ipr_expression = r"IPR[\d]{6}[^0-9]"

    df_genome = df_genome[['Gene_Name', 'GOs', 'EnzymeCodes', 'InterProScan']]
    df_genome.set_index("Gene_Name", inplace = True)
    df_genome['InterProScan'] = df_genome['InterProScan'].str.split("; ").apply(
        lambda x: [y[:len('IPRXXXXXX')] 
                   for y in x 
                   if re.match(ipr_expression, y)])

def ec_extraction(df_genome):
    df_genome['EnzymeCodes'] = df_genome['EnzymeCodes'].str.lower().str.split(";")

    ecs_requests = []
    for ecs in df_genome['EnzymeCodes']:
        if type(ecs) != float:
            for ec in ecs:
                if ec not in ecs_requests:
                    ecs_requests.append(ec)

    return ecs_requests

def r_keggrest_ec(ecs_requests):
    '''
        Uses the R script to extract pathway using Enzyme Code in KEGG.
    '''
    command = 'Rscript'
    path_script = 'pathway_extraction/keggrest_pathway_extraction.R'
    data_name = ["enzyme"]
    cmd = [command, path_script] + ecs_requests + data_name
    try:
        subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as e:
        raise RuntimeError("command '{}' return with error (code {}): {}".format(e.cmd, e.returncode, e.output))

def main():
    if os.path.exists(temporary_directory_database[:-1]) == False :
        os.makedirs(temporary_directory_database)
    name_reference_file = 'Annotation_blast2go_PROT_eH'
    df_genome = pa.read_csv(temporary_directory + name_reference_file + 'GOsTranslatedAndFixed.tsv', sep = "\t")
    ecs_requests = ec_extraction(df_genome)
    r_keggrest_ec(ecs_requests)

    chebi_from_go.go_to_chebi()

    database_mapping_from_gos.main()

    eupathdb_pathway_extraction.main()

    ghost_koala_pathway_extraction.main()

    interpro_pathway_extraction.main()

    panther_pathway_mapping_uniprot.main()

    reactome_pathway_extraction.main(name_reference_file)

    sparql_query_reactome_pathway_name.main()
main()

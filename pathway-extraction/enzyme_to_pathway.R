#!/usr/bin/RenzymeToPathway

library(KEGGREST)

my_args <- commandArgs(trailingOnly = TRUE)

enzyme_codes = my_args

enzyme_pathway_table <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(enzyme_pathway_table) <- c('ecCode', 'Pathway', 'PathwayID')

write.table(enzyme_pathway_table, file <- file.path('temporaryFiles', 'databases', 'enzyme_pathway_kegg.tsv'), append = FALSE, row.names = FALSE, col.names = TRUE, sep ="\t")

for(enzyme_code in enzyme_codes)
{
    pathway_linked_to_enzyme <- keggLink('pathway', enzyme_code)
    if (toString(pathway_linked_to_enzyme) == '')
    {
        enzyme_pathway_table <- matrix(c(enzyme_code , NA, NA), ncol = 3)

        write.table(enzyme_pathway_table, file <- file.path('temporaryFiles', 'databases', 'enzyme_pathway_kegg.tsv'), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t")
    }

    for (pathway_id in pathway_linked_to_enzyme)
    {
        pathway_name <- keggFind('pathway', pathway_id)

        if (grepl(toString(pathway_name), "path:ec") == FALSE)
        {
            enzyme_pathway_table <- matrix(c(enzyme_code, pathway_name, pathway_id), ncol = 3)

            write.table(enzyme_pathway_table, file <- file.path('temporaryFiles', 'databases', 'enzyme_pathway_kegg.tsv'), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t")
        }
    }
}

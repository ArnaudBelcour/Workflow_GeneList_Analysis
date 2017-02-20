#!/usr/bin/RenzymeToPathway

library(KEGGREST)

trim <- function (x) gsub("^\\s+|\\s+$", "", x)

my_args <- commandArgs(trailingOnly = TRUE)

data_codes = my_args[-length(my_args)]

data_name = my_args[length(my_args)]

pathway_table <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(pathway_table) <- c('ecCode', 'Pathway', 'PathwayID')

file_name = paste(data_name, '_pathway_kegg.tsv', sep = "")

write.table(pathway_table, file <- file.path('temporaryFiles', 'databases', file_name), append = FALSE, row.names = FALSE, col.names = TRUE, quote = FALSE, sep ="\t")

for(data_code in data_codes)
{
    pathway_linked_to_enzyme <- keggLink('pathway', data_code)
    if (toString(pathway_linked_to_enzyme) == '')
    {
        pathway_table <- matrix(c(trim(data_code), NA, NA), ncol = 3)

        write.table(pathway_table, file <- file.path('temporaryFiles', 'databases', file_name), append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    }

    for (pathway_id in pathway_linked_to_enzyme)
    {
        pathway_name <- keggFind('pathway', pathway_id)

        if (grepl(toString(pathway_name), "path:ec") == FALSE)
        {
            pathway_table <- matrix(c(trim(data_code), pathway_name, pathway_id), ncol = 3)

            write.table(pathway_table, file <- file.path('temporaryFiles', 'databases', file_name), append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
        }
    }
}

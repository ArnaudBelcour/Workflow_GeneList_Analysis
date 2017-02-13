#!/usr/bin/RenzymeToPathway

library(KEGGREST)

myArgs <- commandArgs(trailingOnly = TRUE)

enzymeCodes = myArgs

enzymePathwayTable <- data.frame(matrix(nrow = 0, ncol = 3))
colnames(enzymePathwayTable) <- c('ecCode', 'Pathway', 'PathwayID')

temporaryFilesDatabases <- '../temporaryFiles/databases'

print(getwd())

write.table(enzymePathwayTable, file <- file.path('temporaryFiles', 'databases', 'enzymeToPathway.tsv'), append = FALSE, row.names = FALSE, col.names = TRUE, sep ="\t")

for(enzymeCode in enzymeCodes)
{
    pathwayLinkedToEnzyme <- keggLink('pathway', enzymeCode)
    if (toString(pathwayLinkedToEnzyme) == '')
    {
        enzymePathwayTable <- matrix(c(enzymeCode, NA, NA), ncol = 3)

        write.table(enzymePathwayTable, file <- file.path('temporaryFiles', 'databases', 'enzymeToPathway.tsv'), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t")
    }

    for (pathwayID in pathwayLinkedToEnzyme)
    {
        pathwayName <- keggFind('pathway', pathwayID)

        if (grepl(toString(pathwayName), "path:ec") == FALSE)
        {
            enzymePathwayTable <- matrix(c(enzymeCode, pathwayName, pathwayID), ncol = 3)

            write.table(enzymePathwayTable, file <- file.path('temporaryFiles', 'databases', 'enzymeToPathway.tsv'), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t")
        }
    }
}

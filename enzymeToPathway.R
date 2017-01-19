library(KEGGREST)

enzymeCodes = c('ec:3.6.1.3', 'ec:3.2.1.55', 'ec:3.2.1.78', 'ec:3.2.1.8', 'ec:2.7.3')

enzymePathwayTable <- data.frame(matrix(ncol=3))
colnames(enzymePathwayTable) <- c('ecCode', 'Pathway', 'PathwayID')
enzymePathwayTable = enzymePathwayTable[-1,]

temporaryFiles = 'temporaryFiles/'

write.table(enzymePathwayTable, file = paste(temporaryFiles, "enzymeToPathway.tsv"), append = FALSE, row.names = FALSE, col.names = TRUE, sep="\t")

for(enzymeCode in enzymeCodes)
{
    pathwayLinkedToEnzyme = keggLink('pathway', enzymeCode)
    if (toString(pathwayLinkedToEnzyme) == '')
    {
        enzymePathwayTable <- matrix(c(enzymeCode, NA, NA), ncol = 3)

        write.table(enzymePathwayTable, file = paste(temporaryFiles, "enzymeToPathway.tsv"), append = T, col.names = FALSE, row.names = FALSE, sep = "\t")
    }

    for (pathwayID in pathwayLinkedToEnzyme)
    {
        pathwayName = keggFind('pathway', pathwayID)

        if (grepl(toString(pathwayName), "path:ec") == FALSE)
        {
            enzymePathwayTable <- matrix(c(enzymeCode, pathwayName, pathwayID), ncol = 3)

            write.table(enzymePathwayTable, file = paste(temporaryFiles, "enzymeToPathway.tsv"), append = TRUE, col.names = FALSE, row.names = FALSE, sep = "\t")
        }
    }
}

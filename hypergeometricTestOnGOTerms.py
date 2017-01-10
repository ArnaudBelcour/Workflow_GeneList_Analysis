import csv
import pandas as pa
import scipy.stats as stats
from ast import literal_eval
import primaryFileManagement

def createGenGOAnalysisFile(fileName):
    table = pa.read_csv(fileName + ".tsv", sep = "\t")
    table = table[['Nom_Gene', 'GOs']]
    table = table.set_index("Nom_Gene")

    csvfile = open("Gene_GO.tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(('Nom_Gene', 'GO'))
    for gene, GOListObject in table.iterrows():
        for GOList in GOListObject:
            for GO in literal_eval(GOList):
                writer.writerow((gene, GO))

    csvfile.close()

def computeHypergeometricTestForeachValue(GOSetNumber, numberOfGeneGenome, GOGenomeNumber, numberOfGeneOfInterest):
    return stats.hypergeom.sf(GOSetNumber - 1, numberOfGeneGenome, GOGenomeNumber, numberOfGeneOfInterest)


def HypergeometricTestOnDataFrame():
    df = pa.read_csv("Gene_GO.tsv", sep = "\t")
    numberOfGeneOfInterest = len(df["Nom_Gene"].unique())

    counts_df = pa.DataFrame(df.groupby("GO").size().rename("Counts"))
    counts_df = counts_df.sort_values(['Counts'], ascending=[False])
    counts_df.reset_index()

    numberOfAnnotations = counts_df['Counts'].sum()

    computeHypergeometricTest = lambda x: computeHypergeometricTestForeachValue(x, 1070, 20, numberOfGeneOfInterest)
    counts_df['HypergeoTest'] = counts_df['Counts'].apply(computeHypergeometricTest)
    print("BonFerroni p-value : " + str((0.05 / numberOfAnnotations)))
    print(counts_df)

def main():
    primaryFileManagement.main()
    createGenGOAnalysisFile("queryResultsGOTranslatedAndFixed")
    HypergeometricTestOnDataFrame()

main()

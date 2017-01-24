import csv
import math
import pandas as pa
import primaryFileManagement
import scipy.stats as stats
import sys
from ast import literal_eval

inputDirectory = "inputFiles/"
temporaryDirectory = 'temporaryFiles/'
outputDirectory = 'outputFiles/'

def createGenGOAnalysisFile(fileName, columnsNames):
    table = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
    table = table[columnsNames]
    table = table.set_index(columnsNames[0])

    csvfile = open(temporaryDirectory + "Gene_GO.tsv", "w", newline = "")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow((columnsNames[0], columnsNames[1]))

    for gene, GOListObject in table.iterrows():
        for GOList in GOListObject:
            for GO in literal_eval(GOList):
                writer.writerow((gene, GO))

    csvfile.close()

def goCountingGeneList(fileName, columnName):
    df = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
    counts_df = pa.DataFrame(df.groupby('GOs').size().rename(columnName))

    numberOfGene = len(df["Gene_Name"].unique())

    counts_df = pa.DataFrame(df.groupby("GOs").size().rename(columnName))
    counts_df = counts_df.sort_values([columnName], ascending=[False])
    counts_df.reset_index()

    numberOfAnnotations = counts_df[columnName].sum()

    return counts_df, numberOfGene

def goCountingGenome(fileName, columnName):
    df = pa.read_csv(inputDirectory + fileName + ".tsv", sep="\t")
    df.columns = ['GOs']
    counts_df_GOGenome = pa.DataFrame(df.groupby("GOs").size().rename(columnName))

    counts_df_GOGenome.reset_index(inplace=True)
    counts_df_GOGenome['GOs'] = counts_df_GOGenome['GOs'].str.replace(":", "_")
    counts_df_GOGenome = counts_df_GOGenome.set_index('GOs')

    return counts_df_GOGenome

def hypergeometricTestOnDataframe(df, numberOfGeneOfInterest, numberOfGenesInGenome, overOrUnderRepresentation, genomeColumns):
    GOtermsWithHyperGeoTestNAN = []

    if overOrUnderRepresentation == "over":
        for GO, row in df.iterrows():
            if math.isnan(df.get_value(GO, genomeColumns)):
                df = df.drop([GO])
            else:
                computeHypergeometricTestOverRepresentation(GO, numberOfGeneOfInterest, row['Counts'], numberOfGenesInGenome, row[genomeColumns], df)
                if math.isnan(df.get_value(GO, 'pValueHypergeometric')):
                    GOtermsWithHyperGeoTestNAN.append(GO)
                    df = df.drop([GO])

    if overOrUnderRepresentation == "under":
        for GO, row in df.iterrows():
            if math.isnan(df.get_value(GO, genomeColumns)):
                df = df.drop([GO])
            else:
                computeHypergeometricTestUnderRepresentation(GO, numberOfGeneOfInterest, row['Counts'], numberOfGenesInGenome, row[genomeColumns], df)
                if math.isnan(df.get_value(GO, 'pValueHypergeometric')):
                    GOtermsWithHyperGeoTestNAN.append(GO)
                    df = df.drop([GO])

    df = df.sort_values("pValueHypergeometric")

    return df

def computeHypergeometricTestOverRepresentation(GO, numberOfGeneOfInterest, GOSetNumber, numberOfGeneGenome, GOGenomeNumber, df):
    pValueHypergeo = stats.hypergeom.sf(GOSetNumber - 1, (numberOfGeneGenome + GOGenomeNumber), GOGenomeNumber, numberOfGeneOfInterest)
    df.set_value(GO, 'pValueHypergeometric', pValueHypergeo)

    return df

def computeHypergeometricTestUnderRepresentation(GO, numberOfGeneOfInterest, GOSetNumber, numberOfGeneGenome, GOGenomeNumber, df):
    pValueHypergeo = stats.hypergeom.cdf(GOSetNumber - 1, (numberOfGeneGenome + GOGenomeNumber), GOGenomeNumber, numberOfGeneOfInterest)
    df.set_value(GO, 'pValueHypergeometric', pValueHypergeo)

    return df

def countingApproximation(df):
    GOtermsWithHyperGeoTestNAN = []

    for GO, row in df.iterrows():
        df.set_value(GO, 'CountsTotal', row['Counts'] + row['CountsGenome'])

    return df

def percentageCalculation(number1, number2):
    percentage = (number1 / number2) * 100

    return percentage

def multipleTestingCorrectionAndOutputWrting(df, numberOfGeneOfInterest, alpha, overOrUnderRepresentation, approximationYesOrNo, yesAnswers):
    df = df.sort_values(['pValueHypergeometric'])

    df = correctionBonferroni(df, numberOfGeneOfInterest)
    df = correctionBenjaminiHochberg(df, numberOfGeneOfInterest)
    df = correctionHolm(df, numberOfGeneOfInterest)

    d_GOLabelToNumber, d_GOLabelWithSynonym = primaryFileManagement.GOLabelNumberDictionnaryCreation(inputDirectory + "queryResults.csv", 'inverse')

    for GO, row in df.iterrows():
        if GO in d_GOLabelToNumber:
            df.set_value(GO, 'GOLabel', d_GOLabelToNumber[GO])

    if approximationYesOrNo in yesAnswers:
        df = df[['Counts', 'CountsTotal', 'PercentageGOList', 'PercentageGOGenome', 'pValueHypergeometric', 'pValueBonferroni', 'pValueHolm', \
                    'pValueBenjaminiHochberg', 'GOLabel']]
    else:
        df = df[['Counts', 'CountsGenome', 'PercentageGOList', 'PercentageGOGenome', 'pValueHypergeometric', 'pValueBonferroni', 'pValueHolm', \
                    'pValueBenjaminiHochberg', 'GOLabel']]

    if overOrUnderRepresentation == 'over':
        df.to_csv(outputDirectory + "pValuesOfGOTermAndLabel_over.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)
    elif overOrUnderRepresentation == 'under':
        df.to_csv(outputDirectory + "pValuesOfGOTermAndLabel_under.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)

    errorRateSidak = errorRateAdjustementBonferroni(alpha, numberOfGeneOfInterest)
    GOSignificativesSidak = selectionGOTermWithAdjustedErrorRate(errorRateSidak, df)
    GOLabelSignificativesSidak = tranlsationGONumberToGOLabel(GOSignificativesSidak, d_GOLabelToNumber)

    errorRateBonferroni = errorRateAdjustementBonferroni(alpha, numberOfGeneOfInterest)
    GOSignificativesBonferroni = selectionGOTermWithAdjustedErrorRate(errorRateBonferroni, df)
    GOLabelSignificativesBonferroni= tranlsationGONumberToGOLabel(GOSignificativesBonferroni, d_GOLabelToNumber)

    GOSignificativesHolm = selectionGOTermWithAdjustedPValue("Holm", alpha, df)
    GOLabelSignificativesHolm = tranlsationGONumberToGOLabel(GOSignificativesHolm, d_GOLabelToNumber)

    GOSignificativesBenjaminiHochberg = selectionGOTermWithAdjustedPValue("BenjaminiHochberg", alpha, df)
    GOLabelSignificativesBenjaminiAndHochberg = tranlsationGONumberToGOLabel(GOSignificativesBenjaminiHochberg, d_GOLabelToNumber)

    if overOrUnderRepresentation == 'over':
        csvfile = open(outputDirectory + "significativesGO_over.tsv", "w", newline = "")
    elif overOrUnderRepresentation == 'under':
        csvfile = open(outputDirectory + "significativesGO_under.tsv", "w", newline = "")

    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['GOSidak', 'GOBonferroni', 'GOHolm', 'GOBenjaminiHochberg'])

    for index in range(len(GOLabelSignificativesBenjaminiAndHochberg)):
        if index in range (len(GOLabelSignificativesSidak)):
            GOLabelSignificativesSidakValue =  GOLabelSignificativesSidak[index]
        else :
            GOLabelSignificativesSidakValue =  'nan'
        if index in range (len(GOLabelSignificativesBonferroni)):
            GOLabelSignificativesBonferroniValue =  GOLabelSignificativesBonferroni[index]
        else :
            GOLabelSignificativesBonferroniValue =  'nan'
        if index in range (len(GOLabelSignificativesHolm)):
            GOLabelSignificativesHolmValue =  GOLabelSignificativesHolm[index]
        else :
            GOLabelSignificativesHolmValue =  'nan'

        writer.writerow([GOLabelSignificativesSidakValue, GOLabelSignificativesBonferroniValue, GOLabelSignificativesHolmValue, \
        GOLabelSignificativesBenjaminiAndHochberg[index]])

    csvfile.close()

def correctionBonferroni(df, numberOfGeneOfInterest):
    pValueCorrectionBonferroni = lambda x: x * numberOfGeneOfInterest
    df['pValueBonferroni'] = df['pValueHypergeometric'].apply(pValueCorrectionBonferroni)

    return df

def correctionBenjaminiHochberg(df, numberOfGeneOfInterest):
    for GO, row in df.iterrows():
        pValueCorrectionBenjaminiHochberg = row['pValueHypergeometric'] * (numberOfGeneOfInterest/(df.index.get_loc(GO)+1))
        df.set_value(GO, 'pValueBenjaminiHochberg', pValueCorrectionBenjaminiHochberg)

    return df

def correctionHolm(df, numberOfGeneOfInterest):
    for GO, row in df.iterrows():
        pValueCorrectionHolm = row['pValueHypergeometric'] * (numberOfGeneOfInterest - df.index.get_loc(GO))
        df.set_value(GO, 'pValueHolm', pValueCorrectionHolm)

    return df

def errorRateAdjustementBonferroni(alpha, numberOfGeneOfInterest):
    errorRateAdjusted = alpha / numberOfGeneOfInterest

    return errorRateAdjusted

def errorRateAdjustementSidak(alpha, numberOfGeneOfInterest):
    errorRateAdjusted = (1- math.pow((1-alpha), (1 / numberOfGeneOfInterest)))

    return errorRateAdjusted

def selectionGOTermWithAdjustedErrorRate(errorRate, df):
    GOSignificatives = []

    for GO, row in df.iterrows():
        if row['pValueHypergeometric'] < errorRate :
            GOSignificatives.append(GO)

    return GOSignificatives

def selectionGOTermWithAdjustedPValue(methodName, alpha, df):
    GOSignificatives = []

    for GO, row in df.iterrows():
        if row['pValue' + methodName] < alpha :
            GOSignificatives.append(GO)

    return GOSignificatives

def tranlsationGONumberToGOLabel(GONumbers, d_GOLabelToNumber):
    GOLabels = []

    for GONumber in GONumbers:
        if GONumber in d_GOLabelToNumber:
            GOLabels.append(d_GOLabelToNumber[GONumber])

    return GOLabels

def enrichmentAnalysis():
    pythonVersion = sys.version_info
    primaryFileManagement.columnGOCleaning()
    createGenGOAnalysisFile("queryResultsGOTranslatedAndFixed", ['Gene_Name', 'GOs'])

    sentenceChoiceNumberGene = "Enter the number of genes in the genome of your organism : "
    numberOfGenesInGenome = int(primaryFileManagement.inputPythonFormat(sentenceChoiceNumberGene, pythonVersion))

    counts_df, numberOfGene = goCountingGeneList("Gene_GO", 'Counts')
    counts_df_GOGenome = goCountingGenome("GOTermsPlasmoGenome", 'CountsGenome')

    dfJoined = counts_df.join(counts_df_GOGenome)

    sentenceChoiceAlpha = "Enter the alpha risk : "
    alpha = float(primaryFileManagement.inputPythonFormat(sentenceChoiceAlpha, pythonVersion))

    yesAnswers = ['yes', 'y', 'oui', 'o']
    sentenceChoiceNumberGene = "Is this an approximation of the genome? "
    yesOrNo = primaryFileManagement.inputPythonFormat(sentenceChoiceNumberGene, pythonVersion)

    for GO, row in dfJoined.iterrows():
        dfJoined.set_value(GO, 'PercentageGOList', percentageCalculation(row['Counts'], numberOfGene))

    if yesOrNo in yesAnswers:
        dfJoinedApproximation = countingApproximation(dfJoined)
        for GO, row in dfJoinedApproximation.iterrows():
            dfJoinedApproximation.set_value(GO, 'PercentageGOGenome', percentageCalculation(row['CountsTotal'], numberOfGenesInGenome))

        dfJoinedOverRepresentation = hypergeometricTestOnDataframe(dfJoinedApproximation, numberOfGene, numberOfGenesInGenome, "over", 'CountsTotal')
        multipleTestingCorrectionAndOutputWrting(dfJoinedOverRepresentation, numberOfGene, alpha, 'over', yesOrNo, yesAnswers)

        dfJoinedUnderRepresentation = hypergeometricTestOnDataframe(dfJoinedApproximation, numberOfGene, numberOfGenesInGenome, "under", 'CountsTotal')
        multipleTestingCorrectionAndOutputWrting(dfJoinedUnderRepresentation, numberOfGene, alpha, 'under', yesOrNo, yesAnswers)

    else:
        for GO, row in dfJoined.iterrows():
            dfJoined.set_value(GO, 'PercentageGOGenome', percentageCalculation(row['CountsGenome'], numberOfGenesInGenome))

        dfJoinedOverRepresentation = hypergeometricTestOnDataframe(dfJoined, numberOfGene, numberOfGenesInGenome, "over", 'CountsGenome')
        multipleTestingCorrectionAndOutputWrting(dfJoinedOverRepresentation, numberOfGene, alpha, 'over', yesOrNo, yesAnswers)

        dfJoinedUnderRepresentation = hypergeometricTestOnDataframe(dfJoined, numberOfGene, numberOfGenesInGenome, "under", 'CountsGenome')
        multipleTestingCorrectionAndOutputWrting(dfJoinedUnderRepresentation, numberOfGene, alpha, 'under', yesOrNo, yesAnswers)

enrichmentAnalysis()

import csv
import pandas as pa
import scipy.stats as stats
from ast import literal_eval
import math
import primaryFileManagement

inputDirectory = "inputFiles/"
temporaryDirectory = 'temporaryFiles/'
outputDirectory = 'outputFiles/'

def createGenGOAnalysisFile(fileName, columnsNames):
    table = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
    table = table[columnsNames]
    table = table.set_index(columnsNames[0])

    csvfile = open(temporaryDirectory + "Gene_GO.tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow((columnsNames[0], columnsNames[1]))

    for gene, GOListObject in table.iterrows():
        for GOList in GOListObject:
            for GO in literal_eval(GOList):
                writer.writerow((gene, GO))

    csvfile.close()

def computeHypergeometricTestForeachValue(GO, numberOfGeneOfInterest, GOSetNumber, numberOfGeneGenome, GOGenomeNumber, dfJoined):
    pValueHypergeo = stats.hypergeom.sf(GOSetNumber - 1, numberOfGeneGenome, GOGenomeNumber, numberOfGeneOfInterest)
    dfJoined.set_value(GO, 'pValueHypergeometric', pValueHypergeo)

def HypergeometricTestOnDataFrame():
    df = pa.read_csv(temporaryDirectory + "Gene_GO.tsv", sep = "\t")
    counts_df = pa.DataFrame(df.groupby('GOs').size().rename("Counts"))

    numberOfGeneOfInterest = len(df["Row.names"].unique())

    counts_df = pa.DataFrame(df.groupby("GOs").size().rename("Counts"))
    counts_df = counts_df.sort_values(['Counts'], ascending=[False])
    counts_df.reset_index()

    numberOfAnnotations = counts_df['Counts'].sum()

    return counts_df

def correctionBonferroni(df):
    pValueCorrectionBonferroni = lambda x: x * numberOfGeneOfInterest
    df['pValueBonferroni'] = df['pValueHypergeometric'].apply(pValueCorrectionBonferroni)
    return df

def correctionBenjaminiHochberg(df):

    for GO, row in df.iterrows():
        pValueCorrectionBenjaminiHochberg = row['pValueHypergeometric'] * (numberOfGeneOfInterest/(df.index.get_loc(GO)+1))
        df.set_value(GO, 'pValueBenjaminiHochberg', pValueCorrectionBenjaminiHochberg)
    return df

def correctionHolm(df):

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

def tranlsationGONumberToGOLabem(GONumbers):
	GOLabels = []
	d_GOLabelToNumber, d_GOLabelWithSynonym = GOLabelNumberDictionnaryCreation("queryResults.csv", 'inverse')

	for GONumber in GONumbers:
			if GONumber in d_GOLabelToNumber:
				GOLabels.append(d_GOLabelToNumber[GONumber])

	return GOLabels

def enrichmentAnalysis():
	primaryFileManagement.main()
	createGenGOAnalysisFile("queryResultsGOTranslatedAndFixed", ['Row.names', 'GOs'])
	counts_df = HypergeometricTestOnDataFrame()

	df = pa.read_csv("Gene_GO.tsv", sep = "\t")
	numberOfGeneOfInterest = len(df["Row.names"].unique())

	df = pa.read_csv("GOTermsPlasmoGenome.tsv", sep="\t")
	df.columns = ['GOs']
	counts_df_GOGenome = pa.DataFrame(df.groupby("GOs").size().rename("CountsGenome"))

	counts_df = counts_df.drop(['-'])
	counts_df_GOGenome.reset_index(inplace=True)
	counts_df_GOGenome['GOs'] = counts_df_GOGenome['GOs'].str.replace(":", "_")
	counts_df_GOGenome = counts_df_GOGenome.set_index('GOs')

	dfJoined = counts_df.join(counts_df_GOGenome)

	GOtermsWithHyperGeoTestNAN = []

	for GO, row in dfJoined.iterrows():
		if math.isnan(dfJoined.get_value(GO, 'CountsGenome')):
			dfJoined = dfJoined.drop([GO])
		else:
			computeHypergeometricTestForeachValue(GO, numberOfGeneOfInterest, row['Counts'], 1070, row['CountsGenome'], dfJoined)
			if math.isnan(dfJoined.get_value(GO, 'pValueHypergeometric')):
				GOtermsWithHyperGeoTestNAN.append(GO)
				dfJoined = dfJoined.drop([GO])

	dfJoined = dfJoined.sort_values("pValueHypergeometric")

	alpha = 0.05

	for GO, row in dfJoined.iterrows():
		dfJoined.set_value(GO, 'CountsTotal', row['Counts'] + row['CountsGenome'])

	for GO, row in dfJoined.iterrows():
		if math.isnan(dfJoined.get_value(GO, 'CountsTotal')):
			dfJoined = dfJoined.drop([GO])
		else:
			computeHypergeometricTestForeachValue(GO, numberOfGeneOfInterest, row['Counts'], 12000, row['CountsTotal'], dfJoined)
			if math.isnan(dfJoined.get_value(GO, 'pValueHypergeometric')):
				GOtermsWithHyperGeoTestNAN.append(GO)
				dfJoined = dfJoined.drop([GO])

	dfJoined.sort_values(['pValueHypergeometric'])
	dfJoined = correctionBonferroni(dfJoined)
	dfJoined = correctionBenjaminiHochberg(dfJoined)
	dfJoined = correctionHolm(dfJoined)

	d_GOLabelToNumber, d_GOLabelWithSynonym = GOLabelNumberDictionnaryCreation("queryResults.csv", 'inverse')
	for GO, row in dfJoined.iterrows():
			if GO in d_GOLabelToNumber:
				dfJoined.set_value(GO, 'GOLabel', d_GOLabelToNumber[GO])
	dfJoined = dfJoined.sort_values('pValueBenjaminiHochberg')
	dfJoined[['Counts', 'CountsGenome', 'CountsTotal', 'pValueBenjaminiHochberg', 'GOLabel']]

	dfJoined.to_csv(outputDirectory + "output.tsv", sep= "\t", index = True, header = True, quoting = csv.QUOTE_NONE)

    errorRateSidak = errorRateAdjustementBonferroni(alpha, numberOfGeneOfInterest)
    GOSignificativesSidak = selectionGOTermWithAdjustedErrorRate(errorRateSidak, dfMerged)
    GOLabelSignificativesSidak = tranlsationGONumberToGOLabem(GOSignificativesSidak)

    errorRateBonferroni = errorRateAdjustementBonferroni(alpha, numberOfGeneOfInterest)
    GOSignificativesBonferroni = selectionGOTermWithAdjustedErrorRate(errorRateBonferroni, dfMerged)
    GOLabelSignificativesBonferroni= tranlsationGONumberToGOLabem(GOSignificativesBonferroni)

    GOSignificativesHolm = selectionGOTermWithAdjustedPValue("Holm", alpha, dfMerged)
    GOLabelSignificativesHolm = tranlsationGONumberToGOLabem(GOSignificativesHolm)

    GOSignificativesBenjaminiHochberg = selectionGOTermWithAdjustedPValue("BenjaminiHochberg", alpha, dfMerged)
    GOLabelSignificativesBenjaminiAndHochberg = tranlsationGONumberToGOLabem(GOSignificativesBenjaminiHochberg)

    csvfile = open(outputDirectory + "significativesGO.tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['GOSidak', 'GOBonferroni', 'GOHolm', 'GOBenjaminiHochberg'])

    for index in range(len(GOLabelSignificativesBenjaminiAndHochberg)):
		try :
			GOlabelSignificativesSidak =  GOLabelSignificativesSidak[index]
		except :
			GOlabelSignificativesSidak =  nan
		try :
			GOLabelSignificativesBonferroni =  GOLabelSignificativesBonferroni[index]
		except :
			GOLabelSignificativesBonferroni =  nan
		try :
			GOLabelSignificativesHolm =  GOLabelSignificativesHolm[index]
		except :
			GOLabelSignificativesHolm =  nan

		writer.writerow([GOLabelSignificativesSidak, GOLabelSignificativesBonferroni, GOLabelSignificativesHolm, GOLabelSignificativesBenjaminiAndHochberg[index]])
	csvfile.close()
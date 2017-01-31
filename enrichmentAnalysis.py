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

class EnrichmentAnalysis():

    def __init__(self, columnName):
        self.objectToAnalyze = columnName
        self.outputColumns = ['Counts', 'CountsGenome', 'Percentage' + self.getObjectToAnalyze() + 'List', 'Percentage' + self.getObjectToAnalyze() + 'Genome', 'pValueHypergeometric', 'pValueBonferroni', 'pValueHolm', 'pValueSGoF', 'pValueBenjaminiHochberg']

    def getObjectToAnalyze(self):
        return self.objectToAnalyze

    def getOutputColumns(self):
        return self.outputColumns

    def setOutputColumns(self, index, value):
        self.outputColumns[index] = value

    def createGeneObjectAnalysisFile(self, fileName, columnsNames):
        df = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
        df = df[columnsNames]
        df = df.set_index(columnsNames[0])

        csvfile = open(temporaryDirectory + "Gene_" + self.getObjectToAnalyze() + ".tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow((columnsNames[0], columnsNames[1]))

        for gene, listOflistOfAnalyzedObjects in df.iterrows():
            for listOfAnalyzedObjects in listOflistOfAnalyzedObjects:
                for analyzedObject in literal_eval(listOfAnalyzedObjects):
                    writer.writerow((gene, analyzedObject))

        csvfile.close()

    def countingGeneList(self, fileName, columnName):
        df = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
        counts_df = pa.DataFrame(df.groupby(self.getObjectToAnalyze()).size().rename(columnName))

        numberOfGene = len(df["Gene_Name"].unique())

        counts_df = pa.DataFrame(df.groupby(self.getObjectToAnalyze()).size().rename(columnName))
        counts_df = counts_df.sort_values([columnName], ascending=[False])
        counts_df.reset_index()

        numberOfAnnotations = counts_df[columnName].sum()

        return counts_df, numberOfGene

    def countingGenome(self, fileName, columnName):
        df = pa.read_csv(inputDirectory + fileName + ".tsv", sep="\t")
        df.columns = [self.getObjectToAnalyze()]
        counts_df_Genome = pa.DataFrame(df.groupby(self.getObjectToAnalyze()).size().rename(columnName))

        counts_df_Genome.reset_index(inplace=True)
        counts_df_Genome[self.getObjectToAnalyze()] = counts_df_Genome[self.getObjectToAnalyze()].str.replace(":", "_")
        counts_df_Genome = counts_df_Genome.set_index(self.getObjectToAnalyze())

        return counts_df_Genome

    def hypergeometricTestOnDataframe(self, df, numberOfGeneOfInterest, numberOfGenesInGenome, overOrUnderRepresentation, genomeColumns):
        analyzedObjectsWithHyperGeoTestNAN = []

        if overOrUnderRepresentation == "over":
            for analyzedObject, row in df.iterrows():
                if math.isnan(df.get_value(analyzedObject, genomeColumns)):
                    df = df.drop([analyzedObject])
                else:
                    self.computeHypergeometricTestOverRepresentation(analyzedObject, numberOfGeneOfInterest, row['Counts'], numberOfGenesInGenome, row[genomeColumns], df)
                    if math.isnan(df.get_value(analyzedObject, 'pValueHypergeometric')):
                        analyzedObjectsWithHyperGeoTestNAN.append(analyzedObject)
                        df = df.drop([analyzedObject])

        if overOrUnderRepresentation == "under":
            for analyzedObject, row in df.iterrows():
                if math.isnan(df.get_value(analyzedObject, genomeColumns)):
                    df = df.drop([analyzedObject])
                else:
                    self.computeHypergeometricTestUnderRepresentation(analyzedObject, numberOfGeneOfInterest, row['Counts'], numberOfGenesInGenome, row[genomeColumns], df)
                    if math.isnan(df.get_value(analyzedObject, 'pValueHypergeometric')):
                        analyzedObjectsWithHyperGeoTestNAN.append(analyzedObject)
                        df = df.drop([analyzedObject])

        df = df.sort_values("pValueHypergeometric")

        return df

    def computeHypergeometricTestOverRepresentation(self, analyzedObject, numberOfGeneOfInterest, numberOfObjectOfInterest, numberOfGeneInGenome, numberOfObjectInGenome, df):
        pValueHypergeo = stats.hypergeom.sf(numberOfObjectOfInterest - 1, numberOfGeneInGenome, numberOfObjectInGenome, numberOfGeneOfInterest)
        df.set_value(analyzedObject, 'pValueHypergeometric', pValueHypergeo)

        return df

    def computeHypergeometricTestUnderRepresentation(self, analyzedObject, numberOfGeneOfInterest, numberOfObjectOfInterest, numberOfGeneInGenome, numberOfObjectInGenome, df):
        pValueHypergeo = stats.hypergeom.cdf(numberOfObjectOfInterest + 1, numberOfGeneInGenome, numberOfObjectInGenome, numberOfGeneOfInterest)
        df.set_value(analyzedObject, 'pValueHypergeometric', pValueHypergeo)

        return df

    def countingApproximation(self, df):
        for analyzedObject, row in df.iterrows():
            df.set_value(analyzedObject, 'CountsTotal', row['Counts'] + row['CountsGenome'])

        return df

    def percentageCalculation(self, number1, number2):
        percentage = (number1 / number2) * 100

        return percentage

    def multipleTestingCorrection(self, df, numberOfGeneOfInterest, alpha):
        df = df.sort_values(['pValueHypergeometric'])

        df = self.correctionBonferroni(df, numberOfGeneOfInterest)
        df = self.correctionBenjaminiHochberg(df, numberOfGeneOfInterest)
        df = self.correctionHolm(df, numberOfGeneOfInterest)
        df = self.correctionSGoF(df, numberOfGeneOfInterest, alpha)

        significativeObjects = {}

        errorRateSidak = self.errorRateAdjustementBonferroni(alpha, df, numberOfGeneOfInterest)
        objectSignificativesSidak = self.selectionObjectWithAdjustedErrorRate(errorRateSidak, df)
        significativeObjects['Sidak'] = objectSignificativesSidak

        errorRateBonferroni = self.errorRateAdjustementBonferroni(alpha, df, numberOfGeneOfInterest)
        objectSignificativesBonferroni = self.selectionObjectWithAdjustedErrorRate(errorRateBonferroni, df)
        significativeObjects['Bonferroni'] = objectSignificativesBonferroni

        objectSignificativesHolm = self.selectionObjectWithAdjustedPValue("Holm", alpha, df)
        significativeObjects['Holm'] = objectSignificativesHolm

        objectSignificativesSGoF = self.selectionObjectWithSGoF("SGoF", df)
        significativeObjects['SGoF'] = objectSignificativesSGoF

        objectSignificativesBenjaminiHochberg = self.selectionObjectWithAdjustedPValue("BenjaminiHochberg", alpha, df)
        significativeObjects['BenjaminiHochberg'] = objectSignificativesBenjaminiHochberg

        return df, significativeObjects

    def writingOutput(self, df, significativeObjects, overOrUnderRepresentation, approximationYesOrNo, yesAnswers):

        if approximationYesOrNo in yesAnswers:
            self.setOutputColumns(1, "CountsTotal")
            df = df[self.getOutputColumns()]
        else:
            df = df[self.getOutputColumns()]

        if overOrUnderRepresentation == 'over':
            df.to_csv(outputDirectory + "pValuesOf" + self.getObjectToAnalyze() + "TermAndLabel_over.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)
        elif overOrUnderRepresentation == 'under':
            df.to_csv(outputDirectory + "pValuesOf" + self.getObjectToAnalyze() + "TermAndLabel_under.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)

        if overOrUnderRepresentation == 'over':
            csvfile = open(outputDirectory + "significatives" + self.getObjectToAnalyze() + "_over.tsv", "w", newline = "")
        elif overOrUnderRepresentation == 'under':
            csvfile = open(outputDirectory + "significatives" + self.getObjectToAnalyze() + "_under.tsv", "w", newline = "")

        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow([self.getObjectToAnalyze() + 'Sidak', self.getObjectToAnalyze() + 'Bonferroni', self.getObjectToAnalyze() + 'Holm', self.getObjectToAnalyze() + 'SGoF', self.getObjectToAnalyze() + 'BenjaminiHochberg'])

        for index in range(len(significativeObjects['BenjaminiHochberg'])):
            if index in range(len(significativeObjects['Sidak'])):
                objectSignificativesSidakValue =  significativeObjects['Sidak'][index]
            else :
                objectSignificativesSidakValue =  'nan'
            if index in range(len(significativeObjects['Bonferroni'])):
                objectSignificativesBonferroniValue =  significativeObjects['Bonferroni'][index]
            else :
                objectSignificativesBonferroniValue =  'nan'
            if index in range(len(significativeObjects['Holm'])):
                objectSignificativesHolmValue =  significativeObjects['Holm'][index]
            else :
                objectSignificativesHolmValue =  'nan'
            if index in range(len(significativeObjects['SGoF'])):
                objectSignificativesSGoFValue =  significativeObjects['SGoF'][index]
            else :
                objectSignificativesSGoFValue =  'nan'

            writer.writerow([objectSignificativesSidakValue, objectSignificativesBonferroniValue, objectSignificativesHolmValue, \
            objectSignificativesSGoFValue, significativeObjects['BenjaminiHochberg'][index]])

        csvfile.close()

    def correctionBonferroni(self, df, numberOfGeneOfInterest):
        pValueCorrectionBonferroni = lambda x: x * len(df.index)
        df['pValueBonferroni'] = df['pValueHypergeometric'].apply(pValueCorrectionBonferroni)

        return df

    def correctionBenjaminiHochberg(self, df, numberOfGeneOfInterest):
        for analyzedObject, row in df.iterrows():
            pValueCorrectionBenjaminiHochberg = row['pValueHypergeometric'] * (len(df.index)/(df.index.get_loc(analyzedObject)+1))
            df.set_value(analyzedObject, 'pValueBenjaminiHochberg', pValueCorrectionBenjaminiHochberg)

        return df

    def correctionHolm(self, df, numberOfGeneOfInterest):
        for analyzedObject, row in df.iterrows():
            pValueCorrectionHolm = row['pValueHypergeometric'] * (len(df.index) - df.index.get_loc(analyzedObject))
            df.set_value(analyzedObject, 'pValueHolm', pValueCorrectionHolm)

        return df

    def correctionSGoF(self, df, numberOfGeneOfInterest, alpha):
        F = len(df.index) * alpha
        df = df.sort_values("pValueHypergeometric")
        R = 0

        for analyzedObject, row in df.iterrows():
            if row['pValueHypergeometric'] < alpha:
                R += 1

        df = df.reset_index()

        objectSignificatives = []
        rowNumber = 0
        while stats.chi2.sf(R, F) < alpha:
            df.set_value(rowNumber, 'pValueSGoF', 'significant')
            objectSignificatives.append(df.iloc[rowNumber][self.getObjectToAnalyze()])
            R = R - 1
            rowNumber = rowNumber + 1

        df = df.set_index(self.getObjectToAnalyze())

        for analyzedObject, row in df.iterrows():
            try:
                if analyzedObject not in objectSignificatives:
                    df.set_value(analyzedObject, 'pValueSGoF', 'nonSignificant')
            except:
                df.set_value(analyzedObject, 'pValueSGoF', 'nonSignificant')

        return df

    def errorRateAdjustementBonferroni(self, alpha, df, numberOfGeneOfInterest):
        errorRateAdjusted = alpha / len(df.index)

        return errorRateAdjusted

    def errorRateAdjustementSidak(self, alpha, df, numberOfGeneOfInterest):
        errorRateAdjusted = (1- math.pow((1-alpha), (1 / len(df.index))))

        return errorRateAdjusted

    def selectionObjectWithAdjustedErrorRate(self, errorRate, df):
        objectSignificatives = []

        for analyzedObject, row in df.iterrows():
            if row['pValueHypergeometric'] < errorRate :
                objectSignificatives.append(analyzedObject)

        return objectSignificatives

    def selectionObjectWithAdjustedPValue(self, methodName, alpha, df):
        objectSignificatives = []

        for analyzedObject, row in df.iterrows():
            if row['pValue' + methodName] < alpha :
                objectSignificatives.append(analyzedObject)

        return objectSignificatives

    def selectionObjectWithSGoF(self, methodName, df):
        objectSignificatives = []

        for analyzedObject, row in df.iterrows():
            if row['pValue' + methodName] == 'significant' :
                objectSignificatives.append(analyzedObject)

        return objectSignificatives

    def enrichmentAnalysis(self):
        pythonVersion = sys.version_info
        primaryFileManagement.columnGOCleaning()
        self.createGeneObjectAnalysisFile("queryResults" + self.getObjectToAnalyze() + "TranslatedAndFixed", ['Gene_Name', self.getObjectToAnalyze()])

        sentenceChoiceNumberGene = "Enter the number of genes in the genome of your organism : "
        numberOfGenesInGenome = int(primaryFileManagement.inputPythonFormat(sentenceChoiceNumberGene, pythonVersion))

        counts_df, numberOfGene = self.countingGeneList("Gene_" + self.getObjectToAnalyze(), 'Counts')
        counts_df_Genome = self.countingGenome("GOTermsPlasmoGenome", 'CountsGenome')

        dfJoined = counts_df.join(counts_df_Genome)

        sentenceChoiceAlpha = "Enter the alpha risk : "
        alpha = float(primaryFileManagement.inputPythonFormat(sentenceChoiceAlpha, pythonVersion))

        yesAnswers = ['yes', 'y', 'oui', 'o']
        sentenceChoiceNumberGene = "Is this an approximation of the genome? "
        yesOrNo = primaryFileManagement.inputPythonFormat(sentenceChoiceNumberGene, pythonVersion)

        for analyzedObject, row in dfJoined.iterrows():
            dfJoined.set_value(analyzedObject, 'Percentage' + self.getObjectToAnalyze() + 'List', self.percentageCalculation(row['Counts'], numberOfGene))

        if yesOrNo in yesAnswers:
            dfJoinedApproximation = self.countingApproximation(dfJoined)
            for analyzedObject, row in dfJoinedApproximation.iterrows():
                dfJoinedApproximation.set_value(analyzedObject, 'Percentage' + self.getObjectToAnalyze() + 'Genome', self.percentageCalculation(row['CountsTotal'], numberOfGenesInGenome))

            dfJoinedOverRepresentation = self.hypergeometricTestOnDataframe(dfJoinedApproximation, numberOfGene, numberOfGenesInGenome, "over", 'CountsTotal')
            dfJoinedOverRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedOverRepresentation, numberOfGene, alpha)
            self.writingOutput(dfJoinedOverRepresentation, significativeObjects, "over", yesOrNo, yesAnswers)

            dfJoinedUnderRepresentation = self.hypergeometricTestOnDataframe(dfJoinedApproximation, numberOfGene, numberOfGenesInGenome, "under", 'CountsTotal')
            dfJoinedUnderRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedUnderRepresentation, numberOfGene, alpha)
            self.writingOutput(dfJoinedUnderRepresentation, significativeObjects, "under", yesOrNo, yesAnswers)

        else:
            for analyzedObject, row in dfJoined.iterrows():
                dfJoined.set_value(analyzedObject, 'Percentage' + self.getObjectToAnalyze() + 'Genome', self.percentageCalculation(row['CountsGenome'], numberOfGenesInGenome))

            dfJoinedOverRepresentation = self.hypergeometricTestOnDataframe(dfJoined, numberOfGene, numberOfGenesInGenome, "over", 'CountsGenome')
            dfJoinedOverRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedOverRepresentation, numberOfGene, alpha)
            self.writingOutput(dfJoinedOverRepresentation, significativeObjects, "over", yesOrNo, yesAnswers)

            dfJoinedUnderRepresentation = self.hypergeometricTestOnDataframe(dfJoined, numberOfGene, numberOfGenesInGenome, "under", 'CountsGenome')
            dfJoinedUnderRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedUnderRepresentation, numberOfGene, alpha)
            self.writingOutput(dfJoinedUnderRepresentation, significativeObjects, "under", yesOrNo, yesAnswers)


class GOEnrichmentAnalysis(EnrichmentAnalysis):

    def __init__(self, columnName):
        EnrichmentAnalysis.__init__(self, columnName)
        self.outputColumns.append("GOLabel")

    def tranlsationGONumberToGOLabel(self, goNumbers, d_GOLabelToNumber):
        goLabels = []

        for goNumber in goNumbers:
            if goNumber in d_GOLabelToNumber:
                goLabels.append(d_GOLabelToNumber[goNumber])

        return goLabels

    def multipleTestingCorrection(self, df, numberOfGeneOfInterest, alpha):
        df = df.sort_values(['pValueHypergeometric'])

        df = self.correctionBonferroni(df, numberOfGeneOfInterest)
        df = self.correctionBenjaminiHochberg(df, numberOfGeneOfInterest)
        df = self.correctionHolm(df, numberOfGeneOfInterest)
        df = self.correctionSGoF(df, numberOfGeneOfInterest, alpha)

        d_GOLabelToNumber, d_GOLabelWithSynonym = primaryFileManagement.GOLabelNumberDictionnaryCreation(inputDirectory + "queryResults.csv", 'inverse')

        significativeObjects = {}

        errorRateSidak = self.errorRateAdjustementBonferroni(alpha, df, numberOfGeneOfInterest)
        objectSignificativesSidak = self.selectionObjectWithAdjustedErrorRate(errorRateSidak, df)
        goLabelSignificativesSidak = self.tranlsationGONumberToGOLabel(objectSignificativesSidak, d_GOLabelToNumber)
        significativeObjects['Sidak'] = goLabelSignificativesSidak

        errorRateBonferroni = self.errorRateAdjustementBonferroni(alpha, df, numberOfGeneOfInterest)
        objectSignificativesBonferroni = self.selectionObjectWithAdjustedErrorRate(errorRateBonferroni, df)
        goLabelSignificativesBonferroni= self.tranlsationGONumberToGOLabel(objectSignificativesBonferroni, d_GOLabelToNumber)
        significativeObjects['Bonferroni'] = goLabelSignificativesBonferroni

        objectSignificativesHolm = self.selectionObjectWithAdjustedPValue("Holm", alpha, df)
        goLabelSignificativesHolm = self.tranlsationGONumberToGOLabel(objectSignificativesHolm, d_GOLabelToNumber)
        significativeObjects['Holm'] = goLabelSignificativesHolm

        objectSignificativesSGoF = self.selectionObjectWithSGoF("SGoF", df)
        goLabelSignificativesSGoF = self.tranlsationGONumberToGOLabel(objectSignificativesSGoF, d_GOLabelToNumber)
        significativeObjects['SGoF'] = goLabelSignificativesSGoF

        objectSignificativesBenjaminiHochberg = self.selectionObjectWithAdjustedPValue("BenjaminiHochberg", alpha, df)
        goLabelSignificativesBenjaminiAndHochberg = self.tranlsationGONumberToGOLabel(objectSignificativesBenjaminiHochberg, d_GOLabelToNumber)
        significativeObjects['BenjaminiHochberg'] = goLabelSignificativesBenjaminiAndHochberg

        for GO, row in df.iterrows():
            if GO in d_GOLabelToNumber:
                df.set_value(GO, 'GOLabel', d_GOLabelToNumber[GO])

        return df, significativeObjects

goEnrichmentAnalysis = GOEnrichmentAnalysis('GOs')
goEnrichmentAnalysis.enrichmentAnalysis()

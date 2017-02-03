#!/usr/bin/env python3

import csv
import math
import pandas as pa
import scipy.stats as stats
import six
import sys
from ast import literal_eval

inputDirectory = "inputFiles/"
temporaryDirectory = 'temporaryFiles/'
outputDirectory = 'outputFiles/'

class EnrichmentAnalysis():

    def __init__(self, columnName):
        self.objectToAnalyze = columnName
        self.outputColumns = ['Counts', 'CountsGenome', 'Percentage' + self.getObjectToAnalyze() + 'List', 'Percentage' + self.getObjectToAnalyze() + 'Genome', 'pValueHypergeometric', 'pValueBonferroni', 'pValueHolm', 'pValueSGoF', 'pValueBenjaminiHochberg']
        self.fileOfInterest = ""
        self.fileOfGenome = ""
        self.numberOfGeneOfInterest = 0
        self.numberOfGeneInGenome = 0
        self.alpha = 0.00

    def getObjectToAnalyze(self):
        return self.objectToAnalyze

    def getOutputColumns(self):
        return self.outputColumns

    def getFileOfInterest(self):
        return self.fileOfInterest

    def getFileOfGenome(self):
        return self.fileOfGenome

    def getNumberOfGeneOfInterest(self):
        return self.numberOfGeneOfInterest

    def getNumberOfGeneInGenome(self):
        return self.numberOfGeneInGenome

    def getAlpha(self):
        return self.alpha

    def setFileOfInterest(self, fileName):
        self.fileOfInterest = fileName

    def setFileOfGenome(self, fileName):
        self.fileOfGenome = fileName

    def setOutputColumns(self, index, value):
        self.outputColumns[index] = value

    def setNumberOfGeneOfInterest(self, value):
        self.numberOfGeneOfInterest = value

    def setNumberOfGeneInGenome(self, value):
        self.numberOfGeneInGenome = value

    def setAlpha(self, value):
        self.alpha = value

    def hypergeometricTestOnDataframe(self, df, overOrUnderRepresentation, genomeColumns):
        analyzedObjectsWithHyperGeoTestNAN = []

        if overOrUnderRepresentation == "over":
            for analyzedObject, row in df.iterrows():
                if math.isnan(df.get_value(analyzedObject, genomeColumns)):
                    df = df.drop([analyzedObject])
                else:
                    self.computeHypergeometricTestOverRepresentation(analyzedObject, row['Counts'], row[genomeColumns], df)
                    if math.isnan(df.get_value(analyzedObject, 'pValueHypergeometric')):
                        analyzedObjectsWithHyperGeoTestNAN.append(analyzedObject)
                        df = df.drop([analyzedObject])

        if overOrUnderRepresentation == "under":
            for analyzedObject, row in df.iterrows():
                if math.isnan(df.get_value(analyzedObject, genomeColumns)):
                    df = df.drop([analyzedObject])
                else:
                    self.computeHypergeometricTestUnderRepresentation(analyzedObject, row['Counts'], row[genomeColumns], df)
                    if math.isnan(df.get_value(analyzedObject, 'pValueHypergeometric')):
                        analyzedObjectsWithHyperGeoTestNAN.append(analyzedObject)
                        df = df.drop([analyzedObject])

        df = df.sort_values("pValueHypergeometric")

        return df

    def computeHypergeometricTestOverRepresentation(self, analyzedObject, numberOfObjectOfInterest, numberOfObjectInGenome, df):
        pValueHypergeo = stats.hypergeom.sf(numberOfObjectOfInterest - 1, self.getNumberOfGeneInGenome(), numberOfObjectInGenome, self.getNumberOfGeneOfInterest())
        df.set_value(analyzedObject, 'pValueHypergeometric', pValueHypergeo)

        return df

    def computeHypergeometricTestUnderRepresentation(self, analyzedObject, numberOfObjectOfInterest, numberOfObjectInGenome, df):
        pValueHypergeo = stats.hypergeom.cdf(numberOfObjectOfInterest + 1, self.getNumberOfGeneInGenome(), numberOfObjectInGenome, self.getNumberOfGeneOfInterest())
        df.set_value(analyzedObject, 'pValueHypergeometric', pValueHypergeo)

        return df

    def countingApproximation(self, df):
        for analyzedObject, row in df.iterrows():
            df.set_value(analyzedObject, 'CountsTotal', row['Counts'] + row['CountsGenome'])

        return df

    def percentageCalculation(self, number1, number2):
        percentage = (number1 / number2) * 100

        return percentage

    def multipleTestingCorrection(self, df):
        df = df.sort_values(['pValueHypergeometric'])

        df = self.correctionBonferroni(df)
        df = self.correctionBenjaminiHochberg(df)
        df = self.correctionHolm(df)
        df = self.correctionSGoF(df)

        significativeObjects = {}

        errorRateSidak = self.errorRateAdjustementBonferroni(df)
        objectSignificativesSidak = self.selectionObjectWithAdjustedErrorRate(errorRateSidak, df)
        significativeObjects['Sidak'] = objectSignificativesSidak

        errorRateBonferroni = self.errorRateAdjustementBonferroni(df)
        objectSignificativesBonferroni = self.selectionObjectWithAdjustedErrorRate(errorRateBonferroni, df)
        significativeObjects['Bonferroni'] = objectSignificativesBonferroni

        objectSignificativesHolm = self.selectionObjectWithAdjustedPValue("Holm", df)
        significativeObjects['Holm'] = objectSignificativesHolm

        objectSignificativesSGoF = self.selectionObjectWithSGoF("SGoF", df)
        significativeObjects['SGoF'] = objectSignificativesSGoF

        objectSignificativesBenjaminiHochberg = self.selectionObjectWithAdjustedPValue("BenjaminiHochberg", df)
        significativeObjects['BenjaminiHochberg'] = objectSignificativesBenjaminiHochberg

        return df, significativeObjects

    def writingOutput(self, df, significativeObjects, overOrUnderRepresentation, approximationYesOrNo, yesAnswers):
        df = df.sort_values(['pValueBenjaminiHochberg'])

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

    def correctionBonferroni(self, df):
        pValueCorrectionBonferroni = lambda x: x * len(df.index)
        df['pValueBonferroni'] = df['pValueHypergeometric'].apply(pValueCorrectionBonferroni)

        return df

    def correctionBenjaminiHochberg(self, df):
        for analyzedObject, row in df.iterrows():
            pValueCorrectionBenjaminiHochberg = row['pValueHypergeometric'] * (len(df.index)/(df.index.get_loc(analyzedObject)+1))
            df.set_value(analyzedObject, 'pValueBenjaminiHochberg', pValueCorrectionBenjaminiHochberg)

        return df

    def correctionHolm(self, df):
        for analyzedObject, row in df.iterrows():
            pValueCorrectionHolm = row['pValueHypergeometric'] * (len(df.index) - df.index.get_loc(analyzedObject))
            df.set_value(analyzedObject, 'pValueHolm', pValueCorrectionHolm)

        return df

    def correctionSGoF(self, df):
        F = len(df.index) * self.getAlpha()
        df = df.sort_values("pValueHypergeometric")
        R = 0

        for analyzedObject, row in df.iterrows():
            if row['pValueHypergeometric'] < self.getAlpha():
                R += 1

        df = df.reset_index()

        objectSignificatives = []
        rowNumber = 0
        while stats.chi2.sf(R, F) < self.getAlpha():
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

    def errorRateAdjustementBonferroni(self, df):
        errorRateAdjusted = self.getAlpha() / len(df.index)

        return errorRateAdjusted

    def errorRateAdjustementSidak(self, df):
        errorRateAdjusted = (1- math.pow((1-self.getAlpha()), (1 / len(df.index))))

        return errorRateAdjusted

    def selectionObjectWithAdjustedErrorRate(self, errorRate, df):
        objectSignificatives = []

        for analyzedObject, row in df.iterrows():
            if row['pValueHypergeometric'] < errorRate :
                objectSignificatives.append(analyzedObject)

        return objectSignificatives

    def selectionObjectWithAdjustedPValue(self, methodName, df):
        objectSignificatives = []

        for analyzedObject, row in df.iterrows():
            if row['pValue' + methodName] < self.getAlpha() :
                objectSignificatives.append(analyzedObject)

        return objectSignificatives

    def selectionObjectWithSGoF(self, methodName, df):
        objectSignificatives = []

        for analyzedObject, row in df.iterrows():
            if row['pValue' + methodName] == 'significant' :
                objectSignificatives.append(analyzedObject)

        return objectSignificatives

    def enrichmentAnalysis(self):

        counts_df = pa.read_csv(temporaryDirectory + self.getFileOfInterest() + ".tsv", sep = "\t")
        counts_df_Genome = pa.read_csv(temporaryDirectory + self.getFileOfGenome() + ".tsv", sep = "\t")

        counts_df = counts_df.set_index("GOs")
        counts_df_Genome = counts_df_Genome.set_index("GOs")

        dfJoined = counts_df.join(counts_df_Genome)
        dfJoined = dfJoined.reset_index()

        yesAnswers = ['yes', 'y', 'oui', 'o']
        sentenceChoiceNumberGene = "Is this an approximation of the genome? "
        yesOrNo = input(sentenceChoiceNumberGene)

        for analyzedObject, row in dfJoined.iterrows():
            dfJoined.set_value(analyzedObject, 'Percentage' + self.getObjectToAnalyze() + 'List', self.percentageCalculation(row['Counts'], self.getNumberOfGeneOfInterest()))

        if yesOrNo in yesAnswers:
            dfJoinedApproximation = self.countingApproximation(dfJoined)
            for analyzedObject, row in dfJoinedApproximation.iterrows():
                dfJoinedApproximation.set_value(analyzedObject, 'Percentage' + self.getObjectToAnalyze() + 'Genome', self.percentageCalculation(row['CountsTotal'], self.getNumberOfGeneInGenome()))

            dfJoinedOverRepresentation = self.hypergeometricTestOnDataframe(dfJoinedApproximation, "over", 'CountsTotal')
            dfJoinedOverRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedOverRepresentation)
            self.writingOutput(dfJoinedOverRepresentation, significativeObjects, "over", yesOrNo, yesAnswers)

            dfJoinedUnderRepresentation = self.hypergeometricTestOnDataframe(dfJoinedApproximation, "under", 'CountsTotal')
            dfJoinedUnderRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedUnderRepresentation)
            self.writingOutput(dfJoinedUnderRepresentation, significativeObjects, "under", yesOrNo, yesAnswers)

        else:
            for analyzedObject, row in dfJoined.iterrows():
                dfJoined.set_value(analyzedObject, 'Percentage' + self.getObjectToAnalyze() + 'Genome', self.percentageCalculation(row['CountsGenome'], self.getNumberOfGeneInGenome()))

            dfJoinedOverRepresentation = self.hypergeometricTestOnDataframe(dfJoined, "over", 'CountsGenome')
            dfJoinedOverRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedOverRepresentation)
            self.writingOutput(dfJoinedOverRepresentation, significativeObjects, "over", yesOrNo, yesAnswers)

            dfJoinedUnderRepresentation = self.hypergeometricTestOnDataframe(dfJoined, "under", 'CountsGenome')
            dfJoinedUnderRepresentation, significativeObjects = self.multipleTestingCorrection(dfJoinedUnderRepresentation)
            self.writingOutput(dfJoinedUnderRepresentation, significativeObjects, "under", yesOrNo, yesAnswers)


class GOEnrichmentAnalysis(EnrichmentAnalysis):

    def __init__(self, columnName, inputDEFileManagement):
        EnrichmentAnalysis.__init__(self, columnName)
        self.outputColumns.append("GOLabel")
        self.inputDEFileInstance = inputDEFileManagement

    def getInputDEFileInstance(self):
        return self.inputDEFileInstance

    def tranlsationGONumberToGOLabel(self, goNumbers, d_GOLabelToNumber):
        goLabels = []

        for goNumber in goNumbers:
            if goNumber in d_GOLabelToNumber:
                goLabels.append(d_GOLabelToNumber[goNumber])

        return goLabels

    def multipleTestingCorrection(self, df):
        df = df.sort_values(['pValueHypergeometric'])

        df = self.correctionBonferroni(df)
        df = self.correctionBenjaminiHochberg(df)
        df = self.correctionHolm(df)
        df = self.correctionSGoF(df)

        d_GOLabelToNumber, d_GOLabelWithSynonym = self.getInputDEFileInstance().GOLabelNumberDictionnaryCreation(inputDirectory + "queryResults.csv", 'inverse')

        significativeObjects = {}

        errorRateSidak = self.errorRateAdjustementBonferroni(df)
        objectSignificativesSidak = self.selectionObjectWithAdjustedErrorRate(errorRateSidak, df)
        goLabelSignificativesSidak = self.tranlsationGONumberToGOLabel(objectSignificativesSidak, d_GOLabelToNumber)
        significativeObjects['Sidak'] = goLabelSignificativesSidak

        errorRateBonferroni = self.errorRateAdjustementBonferroni(df)
        objectSignificativesBonferroni = self.selectionObjectWithAdjustedErrorRate(errorRateBonferroni, df)
        goLabelSignificativesBonferroni= self.tranlsationGONumberToGOLabel(objectSignificativesBonferroni, d_GOLabelToNumber)
        significativeObjects['Bonferroni'] = goLabelSignificativesBonferroni

        objectSignificativesHolm = self.selectionObjectWithAdjustedPValue("Holm", df)
        goLabelSignificativesHolm = self.tranlsationGONumberToGOLabel(objectSignificativesHolm, d_GOLabelToNumber)
        significativeObjects['Holm'] = goLabelSignificativesHolm

        objectSignificativesSGoF = self.selectionObjectWithSGoF("SGoF", df)
        goLabelSignificativesSGoF = self.tranlsationGONumberToGOLabel(objectSignificativesSGoF, d_GOLabelToNumber)
        significativeObjects['SGoF'] = goLabelSignificativesSGoF

        objectSignificativesBenjaminiHochberg = self.selectionObjectWithAdjustedPValue("BenjaminiHochberg", df)
        goLabelSignificativesBenjaminiAndHochberg = self.tranlsationGONumberToGOLabel(objectSignificativesBenjaminiHochberg, d_GOLabelToNumber)
        significativeObjects['BenjaminiHochberg'] = goLabelSignificativesBenjaminiAndHochberg

        for GO, row in df.iterrows():
            if GO in d_GOLabelToNumber:
                df.set_value(GO, 'GOLabel', d_GOLabelToNumber[GO])

        return df, significativeObjects

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

    def getObjectToAnalyze(self):
        return self.objectToAnalyze

    def createGeneObjectAnalysisFile(self, fileName, columnsNames):
        df = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
        df = df[columnsNames]
        df = df.set_index(columnsNames[0])

        csvfile = open(temporaryDirectory + "Gene_" + self.getObjectToAnalyze() + ".tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow((columnsNames[0], columnsNames[1]))

        for gene, listOflistOfBiologicalDatas in df.iterrows():
            for listOfBiologicalDatas in listOflistOfBiologicalDatas:
                for biologicalData in literal_eval(listOfBiologicalDatas):
                    writer.writerow((gene, biologicalData))

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
        GOtermsWithHyperGeoTestNAN = []

        if overOrUnderRepresentation == "over":
            for GO, row in df.iterrows():
                if math.isnan(df.get_value(GO, genomeColumns)):
                    df = df.drop([GO])
                else:
                    self.computeHypergeometricTestOverRepresentation(GO, numberOfGeneOfInterest, row['Counts'], numberOfGenesInGenome, row[genomeColumns], df)
                    if math.isnan(df.get_value(GO, 'pValueHypergeometric')):
                        GOtermsWithHyperGeoTestNAN.append(GO)
                        df = df.drop([GO])

        if overOrUnderRepresentation == "under":
            for GO, row in df.iterrows():
                if math.isnan(df.get_value(GO, genomeColumns)):
                    df = df.drop([GO])
                else:
                    self.computeHypergeometricTestUnderRepresentation(GO, numberOfGeneOfInterest, row['Counts'], numberOfGenesInGenome, row[genomeColumns], df)
                    if math.isnan(df.get_value(GO, 'pValueHypergeometric')):
                        GOtermsWithHyperGeoTestNAN.append(GO)
                        df = df.drop([GO])

        df = df.sort_values("pValueHypergeometric")

        return df

    def computeHypergeometricTestOverRepresentation(self, GO, numberOfGeneOfInterest, GOSetNumber, numberOfGeneInGenome, GOGenomeNumber, df):
        pValueHypergeo = stats.hypergeom.sf(GOSetNumber - 1, numberOfGeneInGenome, GOGenomeNumber, numberOfGeneOfInterest)
        df.set_value(GO, 'pValueHypergeometric', pValueHypergeo)

        return df

    def computeHypergeometricTestUnderRepresentation(self, GO, numberOfGeneOfInterest, GOSetNumber, numberOfGeneInGenome, GOGenomeNumber, df):
        pValueHypergeo = stats.hypergeom.cdf(GOSetNumber + 1, numberOfGeneInGenome, GOGenomeNumber, numberOfGeneOfInterest)
        df.set_value(GO, 'pValueHypergeometric', pValueHypergeo)

        return df

    def countingApproximation(self, df):
        GOtermsWithHyperGeoTestNAN = []

        for GO, row in df.iterrows():
            df.set_value(GO, 'CountsTotal', row['Counts'] + row['CountsGenome'])

        return df

    def percentageCalculation(self, number1, number2):
        percentage = (number1 / number2) * 100

        return percentage

    def multipleTestingCorrectionAndOutputWrting(self, df, numberOfGeneOfInterest, alpha, overOrUnderRepresentation, approximationYesOrNo, yesAnswers):
        df = df.sort_values(['pValueHypergeometric'])

        df = self.correctionBonferroni(df, numberOfGeneOfInterest)
        df = self.correctionBenjaminiHochberg(df, numberOfGeneOfInterest)
        df = self.correctionHolm(df, numberOfGeneOfInterest)
        df = self.correctionSGoF(df, numberOfGeneOfInterest, alpha)

        d_GOLabelToNumber, d_GOLabelWithSynonym = primaryFileManagement.GOLabelNumberDictionnaryCreation(inputDirectory + "queryResults.csv", 'inverse')

        if self.getObjectToAnalyze() == 'GOs':
            for GO, row in df.iterrows():
                if GO in d_GOLabelToNumber:
                    df.set_value(GO, 'GOLabel', d_GOLabelToNumber[GO])

            outputColumns = ['Counts', 'CountsGenome', 'Percentage' + self.getObjectToAnalyze() + 'List', 'Percentage' + self.getObjectToAnalyze() + 'Genome', 'pValueHypergeometric', 'pValueBonferroni', 'pValueHolm', \
                            'pValueSGoF', 'pValueBenjaminiHochberg', 'GOLabel']

            if approximationYesOrNo in yesAnswers:
                outputColumns[1] = 'CountsTotal'
                df = df[outputColumns]
            else:
                df = df[outputColumns]

        else:
            outputColumns = ['Counts', 'CountsGenome', 'Percentage' + self.getObjectToAnalyze() + 'List', 'Percentage' + self.getObjectToAnalyze() + 'Genome', 'pValueHypergeometric', 'pValueBonferroni', 'pValueHolm', \
                            'pValueSGoF', 'pValueBenjaminiHochberg']

            if approximationYesOrNo in yesAnswers:
                outputColumns[1] = 'CountsTotal'
                df = df[outputColumns]
            else:
                df = df[outputColumns]

        if overOrUnderRepresentation == 'over':
            df.to_csv(outputDirectory + "pValuesOf" + self.getObjectToAnalyze() + "TermAndLabel_over.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)
        elif overOrUnderRepresentation == 'under':
            df.to_csv(outputDirectory + "pValuesOf" + self.getObjectToAnalyze() + "TermAndLabel_under.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)

        errorRateSidak = self.errorRateAdjustementBonferroni(alpha, df, numberOfGeneOfInterest)
        objectSignificativesSidak = self.selectionObjectWithAdjustedErrorRate(errorRateSidak, df)

        errorRateBonferroni = self.errorRateAdjustementBonferroni(alpha, df, numberOfGeneOfInterest)
        objectSignificativesBonferroni = self.selectionObjectWithAdjustedErrorRate(errorRateBonferroni, df)


        objectSignificativesHolm = self.selectionObjectWithAdjustedPValue("Holm", alpha, df)
        GOLabelSignificativesHolm = self.tranlsationGONumberToGOLabel(objectSignificativesHolm, d_GOLabelToNumber)

        objectSignificativesBenjaminiHochberg = self.selectionObjectWithAdjustedPValue("BenjaminiHochberg", alpha, df)

        objectSignificativesSGoF = self.selectionObjectWithSGoF("SGoF", df)

        if self.getObjectToAnalyze() == 'GOs':
            GOLabelSignificativesSidak = self.tranlsationGONumberToGOLabel(objectSignificativesSidak, d_GOLabelToNumber)
            GOLabelSignificativesBonferroni= self.tranlsationGONumberToGOLabel(objectSignificativesBonferroni, d_GOLabelToNumber)
            GOLabelSignificativesBenjaminiAndHochberg = self.tranlsationGONumberToGOLabel(objectSignificativesBenjaminiHochberg, d_GOLabelToNumber)
            GOLabelSignificativesSGoF = self.tranlsationGONumberToGOLabel(objectSignificativesSGoF, d_GOLabelToNumber)

        if overOrUnderRepresentation == 'over':
            csvfile = open(outputDirectory + "significatives" + self.getObjectToAnalyze() + "_over.tsv", "w", newline = "")
        elif overOrUnderRepresentation == 'under':
            csvfile = open(outputDirectory + "significatives" + self.getObjectToAnalyze() + "_under.tsv", "w", newline = "")

        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow([self.getObjectToAnalyze() + 'Sidak', self.getObjectToAnalyze() + 'Bonferroni', self.getObjectToAnalyze() + 'Holm', self.getObjectToAnalyze() + 'SGoF', self.getObjectToAnalyze() + 'BenjaminiHochberg'])

        for index in range(len(GOLabelSignificativesBenjaminiAndHochberg)):
            if index in range(len(GOLabelSignificativesSidak)):
                GOLabelSignificativesSidakValue =  GOLabelSignificativesSidak[index]
            else :
                GOLabelSignificativesSidakValue =  'nan'
            if index in range(len(GOLabelSignificativesBonferroni)):
                GOLabelSignificativesBonferroniValue =  GOLabelSignificativesBonferroni[index]
            else :
                GOLabelSignificativesBonferroniValue =  'nan'
            if index in range(len(GOLabelSignificativesHolm)):
                GOLabelSignificativesHolmValue =  GOLabelSignificativesHolm[index]
            else :
                GOLabelSignificativesHolmValue =  'nan'
            if index in range(len(GOLabelSignificativesSGoF)):
                GOLabelSignificativesSGoFValue =  GOLabelSignificativesSGoF[index]
            else :
                GOLabelSignificativesSGoFValue =  'nan'

            writer.writerow([GOLabelSignificativesSidakValue, GOLabelSignificativesBonferroniValue, GOLabelSignificativesHolmValue, \
            GOLabelSignificativesSGoFValue, GOLabelSignificativesBenjaminiAndHochberg[index]])

        csvfile.close()

    def correctionBonferroni(self, df, numberOfGeneOfInterest):
        pValueCorrectionBonferroni = lambda x: x * len(df.index)
        df['pValueBonferroni'] = df['pValueHypergeometric'].apply(pValueCorrectionBonferroni)

        return df

    def correctionBenjaminiHochberg(self, df, numberOfGeneOfInterest):
        for GO, row in df.iterrows():
            pValueCorrectionBenjaminiHochberg = row['pValueHypergeometric'] * (len(df.index)/(df.index.get_loc(GO)+1))
            df.set_value(GO, 'pValueBenjaminiHochberg', pValueCorrectionBenjaminiHochberg)

        return df

    def correctionHolm(self, df, numberOfGeneOfInterest):
        for GO, row in df.iterrows():
            pValueCorrectionHolm = row['pValueHypergeometric'] * (len(df.index) - df.index.get_loc(GO))
            df.set_value(GO, 'pValueHolm', pValueCorrectionHolm)

        return df

    def correctionSGoF(self, df, numberOfGeneOfInterest, alpha):
        F = len(df.index) * alpha
        df = df.sort_values("pValueHypergeometric")
        R = 0

        for GO, row in df.iterrows():
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

        for GO, row in df.iterrows():
            try:
                if GO not in objectSignificatives:
                    df.set_value(GO, 'pValueSGoF', 'nonSignificant')
            except:
                df.set_value(GO, 'pValueSGoF', 'nonSignificant')

        return df

    def errorRateAdjustementBonferroni(self, alpha, df, numberOfGeneOfInterest):
        errorRateAdjusted = alpha / len(df.index)

        return errorRateAdjusted

    def errorRateAdjustementSidak(self, alpha, df, numberOfGeneOfInterest):
        errorRateAdjusted = (1- math.pow((1-alpha), (1 / len(df.index))))

        return errorRateAdjusted

    def selectionObjectWithAdjustedErrorRate(self, errorRate, df):
        objectSignificatives = []

        for GO, row in df.iterrows():
            if row['pValueHypergeometric'] < errorRate :
                objectSignificatives.append(GO)

        return objectSignificatives

    def selectionObjectWithAdjustedPValue(self, methodName, alpha, df):
        objectSignificatives = []

        for GO, row in df.iterrows():
            if row['pValue' + methodName] < alpha :
                objectSignificatives.append(GO)

        return objectSignificatives

    def selectionObjectWithSGoF(self, methodName, df):
        objectSignificatives = []

        for GO, row in df.iterrows():
            if row['pValue' + methodName] == 'significant' :
                objectSignificatives.append(GO)

        return objectSignificatives

    def tranlsationGONumberToGOLabel(self, GONumbers, d_GOLabelToNumber):
        GOLabels = []

        for GONumber in GONumbers:
            if GONumber in d_GOLabelToNumber:
                GOLabels.append(d_GOLabelToNumber[GONumber])

        return GOLabels

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

        for GO, row in dfJoined.iterrows():
            dfJoined.set_value(GO, 'Percentage' + self.getObjectToAnalyze() + 'List', self.percentageCalculation(row['Counts'], numberOfGene))

        if yesOrNo in yesAnswers:
            dfJoinedApproximation = self.countingApproximation(dfJoined)
            for GO, row in dfJoinedApproximation.iterrows():
                dfJoinedApproximation.set_value(GO, 'Percentage' + self.getObjectToAnalyze() + 'Genome', self.percentageCalculation(row['CountsTotal'], numberOfGenesInGenome))

            dfJoinedOverRepresentation = self.hypergeometricTestOnDataframe(dfJoinedApproximation, numberOfGene, numberOfGenesInGenome, "over", 'CountsTotal')
            self.multipleTestingCorrectionAndOutputWrting(dfJoinedOverRepresentation, numberOfGene, alpha, 'over', yesOrNo, yesAnswers)

            dfJoinedUnderRepresentation = self.hypergeometricTestOnDataframe(dfJoinedApproximation, numberOfGene, numberOfGenesInGenome, "under", 'CountsTotal')
            self.multipleTestingCorrectionAndOutputWrting(dfJoinedUnderRepresentation, numberOfGene, alpha, 'under', yesOrNo, yesAnswers)

        else:
            for GO, row in dfJoined.iterrows():
                dfJoined.set_value(GO, 'Percentage' + self.getObjectToAnalyze() + 'Genome', self.percentageCalculation(row['CountsGenome'], numberOfGenesInGenome))

            dfJoinedOverRepresentation = self.hypergeometricTestOnDataframe(dfJoined, numberOfGene, numberOfGenesInGenome, "over", 'CountsGenome')
            self.multipleTestingCorrectionAndOutputWrting(dfJoinedOverRepresentation, numberOfGene, alpha, 'over', yesOrNo, yesAnswers)

            dfJoinedUnderRepresentation = self.hypergeometricTestOnDataframe(dfJoined, numberOfGene, numberOfGenesInGenome, "under", 'CountsGenome')
            self.multipleTestingCorrectionAndOutputWrting(dfJoinedUnderRepresentation, numberOfGene, alpha, 'under', yesOrNo, yesAnswers)

goEnrichmentAnalysis = EnrichmentAnalysis('GOs')
goEnrichmentAnalysis.enrichmentAnalysis()
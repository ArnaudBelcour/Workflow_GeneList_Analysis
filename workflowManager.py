#!/usr/bin/env python3

import os

import primaryFileManagement
from enrichmentAnalysis import GOEnrichmentAnalysis

inputDirectory = "inputFiles/"
temporaryDirectory = 'temporaryFiles/'
temporaryDirectoryDatabase = '../temporaryFiles/databases/'
outputDirectory = 'outputFiles/'

def workflowMainager():
    if os.path.exists(inputDirectory[:-1]) == False :
        os.makedirs(inputDirectory)
        sys.exit("No input data, please put your data fiels in inputFiles directory.")
    if os.path.exists(temporaryDirectory[:-1]) == False :
        os.makedirs(temporaryDirectory)
        os.makedirs(temporaryDirectoryDatabase)
    if os.path.exists(outputDirectory[:-1]) == False :
        os.makedirs(outputDirectory)
 
    if not os.listdir(inputDirectory):
        sys.exit("No input data, please put your data fiels in inputFiles directory.")

    pythonVersion = sys.version_info

    sentenceChoice = "Do you want to download database datas?"
    yesOrNo = primaryFileManagement.inputPythonFormat(sentenceChoice, pythonVersion)
    yesAnswers = ['yes', 'y', 'oui', 'o']

    if yesOrNo in yesAnswers:
        primaryFileManagement.goAncestorsListOfInterest(self.getObjectToAnalyze())

    primaryFileManagement.columnGOCleaning()
    primaryFileManagement.createGeneObjectAnalysisFile("queryResults" + self.getObjectToAnalyze() + "TranslatedAndFixed", ['Gene_Name', self.getObjectToAnalyze()], self.getObjectToAnalyze())

    sentenceChoiceNumberGene = "Enter the number of genes in the genome of your organism : "
    numberOfGenesInGenome = int(primaryFileManagement.inputPythonFormat(sentenceChoiceNumberGene, pythonVersion))

    fileOfInterestName, numberOfGene = primaryFileManagement.countingGeneList("queryResults" + self.getObjectToAnalyze() + "TranslatedAndFixed", 'Counts', self.getObjectToAnalyze())
    fileOfGenomeName = primaryFileManagement.countingGenome("test_genomeGO", 'CountsGenome', self.getObjectToAnalyze())

    sentenceChoiceAlpha = "Enter the alpha risk : "
    alpha = float(primaryFileManagement.inputPythonFormat(sentenceChoiceAlpha, pythonVersion))

    goEnrichmentAnalysis = GOEnrichmentAnalysis('GOs')
    goEnrichmentAnalysis.setFileOfInterest(fileOfInterestName)
    goEnrichmentAnalysis.setFileOfGenome(fileOfGenomeName)
    goEnrichmentAnalysis.setNumberOfGeneOfInterest(numberOfGene)
    goEnrichmentAnalysis.setNumberOfGeneInGenome(numberOfGenesInGenome)
    goEnrichmentAnalysis.setAlpha(alpha)
    goEnrichmentAnalysis.enrichmentAnalysis()

workflowMainager()
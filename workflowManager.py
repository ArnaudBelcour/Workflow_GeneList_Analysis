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

    goEnrichmentAnalysis = GOEnrichmentAnalysis('GOs')
    objectToAnalyze = goEnrichmentAnalysis.getObjectToAnalyze()
    pythonVersion = goEnrichmentAnalysis.getPythonVersion()

    sentenceChoice = "Do you want to download database datas? "
    yesOrNo = primaryFileManagement.inputPythonFormat(sentenceChoice, pythonVersion)
    yesAnswers = ['yes', 'y', 'oui', 'o']

    if yesOrNo in yesAnswers:
        primaryFileManagement.goAncestorsListOfInterest(objectToAnalyze)

    primaryFileManagement.columnGOCleaning(goEnrichmentAnalysis)
    primaryFileManagement.createGeneObjectAnalysisFile("queryResults" + objectToAnalyze + "TranslatedAndFixed", ['Gene_Name', objectToAnalyze], objectToAnalyze)

    sentenceChoiceNumberGene = "Enter the number of genes in the genome of your organism : "
    numberOfGenesInGenome = int(primaryFileManagement.inputPythonFormat(sentenceChoiceNumberGene, pythonVersion))

    fileOfInterestName, numberOfGene = primaryFileManagement.countingGeneList("queryResults" + objectToAnalyze + "TranslatedAndFixed", 'Counts', objectToAnalyze)
    fileOfGenomeName = primaryFileManagement.countingGenome("test_genomeGO", 'CountsGenome', objectToAnalyze)

    sentenceChoiceAlpha = "Enter the alpha risk : "
    alpha = float(primaryFileManagement.inputPythonFormat(sentenceChoiceAlpha, pythonVersion))

    goEnrichmentAnalysis.setFileOfInterest(fileOfInterestName)
    goEnrichmentAnalysis.setFileOfGenome(fileOfGenomeName)
    goEnrichmentAnalysis.setNumberOfGeneOfInterest(numberOfGene)
    goEnrichmentAnalysis.setNumberOfGeneInGenome(numberOfGenesInGenome)
    goEnrichmentAnalysis.setAlpha(alpha)
    goEnrichmentAnalysis.enrichmentAnalysis()

workflowMainager()

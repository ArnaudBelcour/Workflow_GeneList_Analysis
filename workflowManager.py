#!/usr/bin/env python3

import os
import six

from enrichmentAnalysis import GOEnrichmentAnalysis
from fileManagement import FileManagement

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

    sentenceChoice = "Do you want to download database datas? "
    yesOrNo = input(sentenceChoice)
    yesAnswers = ['yes', 'y', 'oui', 'o']

    if yesOrNo in yesAnswers:
        inputDEFileManagement.goAncestorsListOfInterest(objectToAnalyze)

    sentenceChoice = "Write the name of your input file containing differentially expressed gene : "
    nameDEInputFile = input(sentenceChoice)
    inputDEFileManagement = FileManagement(nameDEInputFile)
    fileDENameInput = inputDEFileManagement.getFileName()

    goEnrichmentAnalysis = GOEnrichmentAnalysis('GOs', inputDEFileManagement)
    objectToAnalyze = goEnrichmentAnalysis.getObjectToAnalyze()

    inputDEFileManagement.columnGOCleaning(goEnrichmentAnalysis)
    inputDEFileManagement.createGeneObjectAnalysisFile(fileDENameInput + objectToAnalyze + "TranslatedAndFixed.tsv", ['Gene_Name', objectToAnalyze], objectToAnalyze)
    inputDEFileManagement.goAncestorsListOfInterest(objectToAnalyze)

    sentenceChoice = "Write the name of your input file containing genome : "
    nameGenomeInputFile = input(sentenceChoice)
    inputGenomeFileManagement = FileManagement(nameGenomeInputFile)
    fileGenomeNameInput = inputGenomeFileManagement.getFileName()

    inputGenomeFileManagement.columnGOCleaning(goEnrichmentAnalysis)
    inputGenomeFileManagement.createGeneObjectAnalysisFile(fileGenomeNameInput + objectToAnalyze + "TranslatedAndFixed.tsv", ['Gene_Name', objectToAnalyze], objectToAnalyze)
    inputGenomeFileManagement.goAncestorsListOfInterest(objectToAnalyze)

    fileOfInterestName, numberOfGene = inputDEFileManagement.countingGeneList(fileDENameInput, 'Counts', objectToAnalyze)
    fileOfGenomeName = inputGenomeFileManagement.countingGenome(fileGenomeNameInput, 'CountsGenome', objectToAnalyze)

    sentenceChoiceNumberGene = "Enter the number of genes in the genome of your organism : "
    numberOfGenesInGenome = int(input(sentenceChoiceNumberGene))

    sentenceChoiceAlpha = "Enter the alpha risk : "
    alpha = float(input(sentenceChoiceAlpha))

    goEnrichmentAnalysis.setFileOfInterest(fileOfInterestName)
    goEnrichmentAnalysis.setFileOfGenome(fileOfGenomeName)
    goEnrichmentAnalysis.setNumberOfGeneOfInterest(numberOfGene)
    goEnrichmentAnalysis.setNumberOfGeneInGenome(numberOfGenesInGenome)
    goEnrichmentAnalysis.setAlpha(alpha)
    goEnrichmentAnalysis.enrichmentAnalysis()

workflowMainager()

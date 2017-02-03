#!/usr/bin/env python3

import csv
import math
import os
import pandas as pa
import six
from ast import literal_eval
from collections import defaultdict

import goTermExtractionUniprot

inputDirectory = "inputFiles/"
temporaryDirectory = 'temporaryFiles/'

class FileManagement():

    def __init__(self, nameOftheFile):
        self.fileName = os.path.splitext(nameOftheFile)[0]
        self.fileExtension = os.path.splitext(nameOftheFile)[1]

    def getFileName(self):
        return self.fileName

    def getFileExtension(self):
        return self.fileExtension

    def setFileExtension(self, extension):
        self.fileExtension = extension

    def GOLabelNumberDictionnaryCreation(self, fileName, specification):
        d_GOLabelToNumber = {}

        with open (fileName, 'r') as file:
            queryResultsFile = file.read()
            queryResultsFile = queryResultsFile.replace(" ,\n", "\n")
            queryResultsFile = queryResultsFile.replace(" , ", "\t")
            queryResultsModified = open(temporaryDirectory + "queryResults.tsv", "w")
            queryResultsModified.write(queryResultsFile)
            queryResultsModified.close()

        queryResultsDataframe = pa.read_csv(temporaryDirectory + "queryResults.tsv", sep = "\t")

        queryResultsDataframe.columns = [["subject", "label", "NarrowSynonym", "BroadSynonym", "RelatedSynonym"]]

        quoteDeletion = lambda x: x.replace('"', '')
        queryResultsDataframe["subject"] = queryResultsDataframe["subject"].apply(quoteDeletion)

        goIsolation = lambda x: x[32:]
        queryResultsDataframe["subject"] = queryResultsDataframe["subject"].apply(goIsolation)

        spaceDeletion = lambda x: x.replace(" ", "")
        queryResultsDataframe["label"] = queryResultsDataframe["label"].apply(spaceDeletion)

        spaceDeletion = lambda x: str(x).replace(" ", "")
        queryResultsDataframe["NarrowSynonym"] = queryResultsDataframe["NarrowSynonym"].apply(spaceDeletion)

        spaceDeletion = lambda x: str(x).replace(" ", "")
        queryResultsDataframe["BroadSynonym"] = queryResultsDataframe["BroadSynonym"].apply(spaceDeletion)

        spaceDeletion = lambda x: str(x).replace(" ", "")
        queryResultsDataframe["RelatedSynonym"] = queryResultsDataframe["RelatedSynonym"].apply(spaceDeletion)

        if specification == "inverse":
            d_GOLabelToNumber = dict(zip(queryResultsDataframe['subject'], queryResultsDataframe['label']))

        else:
            d_GOLabelToNumber = dict(zip(queryResultsDataframe['label'], queryResultsDataframe['subject']))

        d_GOLabelWithSynonym = {}

        queryResultsDataframe = queryResultsDataframe.set_index(queryResultsDataframe['subject'])

        for index, row in queryResultsDataframe.iterrows():
            if row['NarrowSynonym'] not in d_GOLabelWithSynonym and row['BroadSynonym'] not in d_GOLabelWithSynonym\
            and row['RelatedSynonym'] not in d_GOLabelWithSynonym:
                if row['NarrowSynonym'] != 'nan':
                    d_GOLabelWithSynonym[row['NarrowSynonym']] = index
                if row['BroadSynonym'] != 'nan':
                    d_GOLabelWithSynonym[row['BroadSynonym']] = index
                if row['RelatedSynonym'] != 'nan':
                    d_GOLabelWithSynonym[row['RelatedSynonym']] = index

        keysTodelete = []

        for key, values in d_GOLabelWithSynonym.items():
            if values == []:
                keysTodelete.append(key)

        for key in keysTodelete:
            del d_GOLabelWithSynonym[key]

        return d_GOLabelToNumber, d_GOLabelWithSynonym

    def translateGOTerm(self, listOfGOLabel, d_GOLabelToNumber):
        listOfGONumber = []

        if type(listOfGOLabel) is float:
            if math.isnan(listOfGOLabel):
                listOfGOLabel = [str(listOfGOLabel)]

        listOfGOLabel = [GO.replace(" ", "") for GO in listOfGOLabel]

        for GOLabel in listOfGOLabel:
            if GOLabel in d_GOLabelToNumber:
                GONumber = d_GOLabelToNumber[GOLabel]
                listOfGONumber.append(GONumber)
            else:
                listOfGONumber.append(GOLabel)

        return listOfGONumber

    def fixObsoleteGOTerm(self, listOfGOLabelAndNumber, d_GOLabelToNumber, d_GOLabelWithSynonym):
        listOfGONumber = []

        listOfGOLabelAndNumber = [GO.replace(" ", "") for GO in listOfGOLabelAndNumber]

        for GOLabelAndNumber in listOfGOLabelAndNumber:
            if GOLabelAndNumber == "-":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] == "GO":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] != 'GO' and GOLabelAndNumber[0:2] != "-":
                GOLabelObsolete = "obsolete" + GOLabelAndNumber
                if GOLabelObsolete in d_GOLabelToNumber:
                    GONumber = d_GOLabelToNumber[GOLabelObsolete]
                    listOfGONumber.append(GONumber)
                elif GOLabelObsolete in d_GOLabelWithSynonym:
                    GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                    listOfGONumber.append(GONumber)
                else:
                    listOfGONumber.append(GOLabelAndNumber)

        return listOfGONumber

    def fixWrongNorLTerm(self, listOfGOLabelAndNumber, d_GOLabelToNumber, d_GOLabelWithSynonym):
        listOfGONumber = []

        listOfGOLabelAndNumber = [GO.replace(" ", "") for GO in listOfGOLabelAndNumber]

        for GOLabelAndNumber in listOfGOLabelAndNumber:
            if GOLabelAndNumber == "-":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] == "GO":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] != 'GO' and GOLabelAndNumber[0:2] != "-":
                if "N" in GOLabelAndNumber :
                    GOLabelNWrong = GOLabelAndNumber.replace("N", "L")
                    if GOLabelNWrong in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelNWrong]
                        listOfGONumber.append(GONumber)
                    elif GOLabelNWrong in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    else:
                        listOfGONumber.append(GOLabelAndNumber)

                elif "L" in GOLabelAndNumber :
                    GOLabelLWrong = GOLabelAndNumber.replace("L", "N")
                    if GOLabelLWrong in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelLWrong]
                        listOfGONumber.append(GONumber)
                    elif GOLabelLWrong in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)

                    else:
                        listOfGONumber.append(GOLabelAndNumber)
                else:
                    listOfGONumber.append(GOLabelAndNumber)

        return listOfGONumber

    def fixProblemsWithSynonym(self, listOfGOLabelAndNumber, d_GOLabelToNumber, d_GOLabelWithSynonym):
        listOfGONumber = []

        listOfGOLabelAndNumber = [GO.replace(" ", "") for GO in listOfGOLabelAndNumber]

        for GOLabelAndNumber in listOfGOLabelAndNumber:
            if GOLabelAndNumber == "-":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] == "GO":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] != 'GO' and GOLabelAndNumber[0:2] != "-":
                if GOLabelAndNumber in d_GOLabelWithSynonym :
                        GONumber = d_GOLabelWithSynonym[GOLabelAndNumber]
                        listOfGONumber.append(GONumber)
                else:
                    listOfGONumber.append(GOLabelAndNumber)

        return listOfGONumber

    def fixTermsIssue(self, listOfGOLabelAndNumber, d_GOLabelToNumber, d_GOLabelWithSynonym):
        listOfGONumber = []

        listOfGOLabelAndNumber = [GO.replace(" ", "") for GO in listOfGOLabelAndNumber]

        for GOLabelAndNumber in listOfGOLabelAndNumber:
            if GOLabelAndNumber == "-":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] == "GO":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] != 'GO' and GOLabelAndNumber[0:2] != "-":
                if "hydrogen" in GOLabelAndNumber :
                    GOLabelFixed = GOLabelAndNumber.replace("hydrogen", "proton")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "al" in GOLabelAndNumber :
                    GOLabelFixed = GOLabelAndNumber.replace("al", "")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "dependent" in GOLabelAndNumber :
                    GOLabelFixed = GOLabelAndNumber.replace("dependent", "templated")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "sequence-specific" in GOLabelAndNumber :
                    GOLabelFixed = GOLabelAndNumber.replace("sequence-specific", "") + ",sequence-specific"
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "homocysteineS-methyltransferaseactivity" in GOLabelAndNumber:
                    GOLabelFixed = "S-adenosylmethionine-" + GOLabelAndNumber
                    if GOLabelFixed in d_GOLabelToNumber[GOLabelFixed]:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "Cul4-RINGubiquitinligasecomplex" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber[:9] + "E3" + GOLabelAndNumber[9:]
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "organ" in GOLabelAndNumber:
                    GOLabelFixed = "animal" + GOLabelAndNumber
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "mitosis" in GOLabelAndNumber:
                    GOLabelFixed = "mitoticnucleardivision"
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "stimulus" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber[:-len("stimulus")]
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "tailtipmorphogenesis" in GOLabelAndNumber:
                    GOLabelFixed = "nematodemale" + GOLabelAndNumber
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "carboxylesteraseactivity" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber[0:len('carboxyl')] + 'icesterhydrolaseactivity'
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "homophiliccelladhesion" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber + 'viaplasmamembraneadhesionmolecules'
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "LSUrRNAbinding" in GOLabelAndNumber:
                    GOLabelFixed = 'largeribosomalsubunit' + GOLabelAndNumber[0:-len('LSU')]
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "threonylcarbamoyladenosine" in GOLabelAndNumber:
                    GOLabelFixed = 'cyclic' + GOLabelAndNumber
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "intraflagellartransportparticleB" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("flagellar", "ciliary")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "intraflagellartransportparticleA" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("flagellar", "ciliary")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "extracellularvesicularexosome" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("vesicular", "")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "Mphaseofmitoticcellcycle" in GOLabelAndNumber:
                    GOLabelFixed = 'mitoticMphase'
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "ATADPantiporteractivity" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber[0:len('AT')] + 'P:' + GOLabelAndNumber[len('AT'):]
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "ribonucleaseHactivity" in GOLabelAndNumber:
                    GOLabelFixed = "RNA-DNAhybridribonucleaseactivity"
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "methylatedhistoneresiduebinding" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("residue", "")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "ciliumaxonemeassembly" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("axoneme", "")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "telomerictemplateRNAreversetranscriptaseactivity" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("ictemplate", "ase")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "ciliumaxoneme" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("cilium", "")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "autophagicvacuolemembrane" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("icvacuole", "osome")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "autophagicvacuoleassembly" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("icvacuole", "osome")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "methylenetetrahydrofolatereductase(NADPH)activity" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber.replace("P", "(P)")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                elif "cytoplasmictransport" in GOLabelAndNumber:
                    GOLabelFixed = GOLabelAndNumber + ",nursecelltooocyte"
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                else:
                    listOfGONumber.append(GOLabelAndNumber)

        return listOfGONumber

    def fixDashInExcess(self, listOfGOLabelAndNumber, d_GOLabelToNumber, d_GOLabelWithSynonym):
        listOfGONumber = []

        listOfGOLabelAndNumber = [GO.replace(" ", "") for GO in listOfGOLabelAndNumber]

        for GOLabelAndNumber in listOfGOLabelAndNumber:
            if GOLabelAndNumber == "-":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] == "GO":
                listOfGONumber.append(GOLabelAndNumber)
            elif GOLabelAndNumber[0:2] != 'GO' and GOLabelAndNumber[0:2] != "-":
                if "-" in GOLabelAndNumber :
                    GOLabelFixed = GOLabelAndNumber.replace("-", "")
                    if GOLabelFixed in d_GOLabelToNumber:
                        GONumber = d_GOLabelToNumber[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                    elif GOLabelFixed in d_GOLabelWithSynonym:
                        GONumber = d_GOLabelWithSynonym[GOLabelFixed]
                        listOfGONumber.append(GONumber)
                else:
                    listOfGONumber.append(GOLabelAndNumber)

        return listOfGONumber

    def cleaningValue(self, dataframe, value):
        valueDataframe = dataframe[dataframe.GOs.str.match(value) == True]

        dataframe = dataframe.set_index("Gene_Name")
        for index in valueDataframe['Gene_Name'].tolist():
            dataframe = dataframe.drop(index)
        dataframe = dataframe.reset_index()

        return dataframe

    def cleaningNanValue(self, dataframe, column):
        for index, row in dataframe.iterrows():
            if type(row[column]) is float:
                if math.isnan(dataframe.get_value(index, column)):
                    dataframe = dataframe.drop([index])

        return dataframe

    def rewritingFile(self, newtable, fileName):
        newtable.to_csv(temporaryDirectory + fileName, "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

    def findColumnOfInterest(self, df):
        columns = []

        for column in df:
            columns.append(column)

        columnsToDeletes = []

        for column in columns:
            for columnValue in df[column]:
                if "no GO terms" in str(columnValue):
                    if column not in columnsToDeletes:
                        columnsToDeletes.append(column)
                if "[EC:" in str(columnValue):
                    if column not in columnsToDeletes:
                        columnsToDeletes.append(column)

        for columnToDelete in columnsToDeletes:
            columns.remove(columnToDelete)

        goColumn = ""
        ecColumn = ""
        iprColumn = ""

        for column in columns:
            for columnValue in df[column]:
                if type(columnValue) is not str:
                    continue
                if ("P:" in columnValue or "F:" in columnValue) and "GO:" not in columnValue:
                    if goColumn == "":
                        goColumn = column
                    else:
                        break
                if "EC:" in columnValue:
                    if ecColumn == "":
                        ecColumn = column
                    else:
                        break
                if "IPR" in columnValue:
                    if iprColumn == "":
                        iprColumn = column
                    else:
                        break

        return goColumn, ecColumn, iprColumn

    def columnGOCleaning(self, instanceEnrichmentAnalysis):

        nameInputFile = self.getFileName()
        extensionInputFile = self.getFileExtension()

        if extensionInputFile == '.xls':
            resultsDataframe = pa.read_excel(inputDirectory + nameInputFile + extensionInputFile, sep = None, na_values = "")
        else:
            resultsDataframe = pa.read_csv(inputDirectory + nameInputFile + extensionInputFile, sep = None, engine = "python", na_values = "")

        sentenceChoice = "Is the first columns of your file, the column containing gene name? "
        yesOrNo = input(sentenceChoice)

        yesAnswers = ['yes', 'y', 'oui', 'o']
        if yesOrNo.lower() in yesAnswers :
            nameGeneColumn = resultsDataframe.columns[0]
        else :
            sentenceChoice = "Write the name of the column containing the gene names : "
            nameGeneColumn = input(sentenceChoice)

        goColumn, ecColumn, iprColumn = self.findColumnOfInterest(resultsDataframe)
        resultsDataframe = resultsDataframe[[nameGeneColumn, goColumn, ecColumn, iprColumn]]
        resultsDataframe.columns = [['Gene_Name', 'GOs', 'EnzymeCodes', 'InterProScan']]

        resultsDataframe['GOs'] = resultsDataframe['GOs'].str.replace("C:", "")
        resultsDataframe['GOs'] = resultsDataframe['GOs'].str.replace("P:", "")
        resultsDataframe['GOs'] = resultsDataframe['GOs'].str.replace("F:", "")
        resultsDataframe['GOs'] = resultsDataframe['GOs'].str.replace(":", "_")
        resultsDataframe['GOs'] = resultsDataframe['GOs'].str.split(";")

        d_GOLabelToNumber, d_GOLabelWithSynonym = self.GOLabelNumberDictionnaryCreation(inputDirectory + "queryResults.csv", 'normal')

        resultsDataframe = self.cleaningValue(resultsDataframe, '-')
        resultsDataframe = self.cleaningNanValue(resultsDataframe, 'GOs')

        translation = lambda x: self.translateGOTerm(x, d_GOLabelToNumber)
        resultsDataframe['GOs'] = resultsDataframe['GOs'].apply(translation)

        correctionProblemsSynonym = lambda x : self.fixProblemsWithSynonym(x,  d_GOLabelToNumber, d_GOLabelWithSynonym)
        resultsDataframe['GOs'] = resultsDataframe['GOs'].apply(correctionProblemsSynonym)

        correctionObsolete = lambda x: self.fixObsoleteGOTerm(x, d_GOLabelToNumber, d_GOLabelWithSynonym)
        resultsDataframe['GOs'] = resultsDataframe['GOs'].apply(correctionObsolete)

        correctionNorL = lambda x: self.fixWrongNorLTerm(x, d_GOLabelToNumber, d_GOLabelWithSynonym)
        resultsDataframe['GOs'] = resultsDataframe['GOs'].apply(correctionNorL)

        correctionTermsIssue = lambda x : self.fixTermsIssue(x, d_GOLabelToNumber, d_GOLabelWithSynonym)
        resultsDataframe['GOs'] = resultsDataframe['GOs'].apply(correctionTermsIssue)

        correctionDashIssue = lambda x : self.fixDashInExcess(x, d_GOLabelToNumber, d_GOLabelWithSynonym)
        resultsDataframe['GOs'] = resultsDataframe['GOs'].apply(correctionDashIssue)

        self.rewritingFile(resultsDataframe, nameInputFile + "GOsTranslatedAndFixed.tsv")

    def createGeneObjectAnalysisFile(self, fileName, columnsNames, columnAnalyzedObject):
        df = pa.read_csv(temporaryDirectory + fileName, sep = "\t")
        df = df[columnsNames]
        df = df.set_index(columnsNames[0])

        csvfile = open(temporaryDirectory + "Gene_" + fileName + columnAnalyzedObject + ".tsv", "w", newline = "")
        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow((columnsNames[0], columnsNames[1]))

        for gene, listOflistOfAnalyzedObjects in df.iterrows():
            for listOfAnalyzedObjects in listOflistOfAnalyzedObjects:
                for analyzedObject in literal_eval(listOfAnalyzedObjects):
                    writer.writerow((gene, analyzedObject))

        csvfile.close()

    def goAncestorsListOfInterest(self, columnAnalyzedObject):
        nameInputFile = self.getFileName()

        df = pa.read_csv(temporaryDirectory + nameInputFile + "GOsTranslatedAndFixed.tsv", '\t')

        for index, row in df.iterrows():
            row[columnAnalyzedObject] = goTermExtractionUniprot.unionGOAndTheirAncestors(literal_eval(row[columnAnalyzedObject]))

        df.to_csv(temporaryDirectory + nameInputFile + "GOsTranslatedAndFixed.tsv", '\t')

    def countingGeneList(self, fileName, columnName, columnAnalyzedObject):
        analyzedObjects = []
        df = pa.read_csv(temporaryDirectory + "Gene_" + fileName + columnAnalyzedObject + ".tsv", sep = "\t")
        for index, row in df.iterrows():
            for analyzedObject in literal_eval(row[columnAnalyzedObject]):
                analyzedObjects.append(analyzedObject)

        counts_df = pa.DataFrame(analyzedObjects)
        counts_df.columns = [columnAnalyzedObject]
        counts_df = counts_df.groupby(columnAnalyzedObject).size().rename(columnName)
        counts_df = counts_df.to_frame()

        numberOfGene = len(df["Gene_Name"].unique())

        counts_df = counts_df.reset_index()
        self.rewritingFile(counts_df, "countingObjectsInInterest.tsv")

        return "countingObjectsInInterest", numberOfGene

    def countingGenome(self, fileName, columnName, columnAnalyzedObject):
        analyzedObjects = []
        df = pa.read_csv(temporaryDirectory + "Gene_" + fileName + columnAnalyzedObject + ".tsv", sep="\t")

        for index, row in df.iterrows():
            for analyzedObject in literal_eval(row[columnAnalyzedObject]):
                analyzedObjects.append(analyzedObject)

        counts_df_Genome = pa.DataFrame(analyzedObjects)
        counts_df_Genome.columns = [columnAnalyzedObject]
        counts_df_Genome = counts_df_Genome.groupby(columnAnalyzedObject).size().rename(columnName)
        counts_df_Genome = counts_df_Genome.to_frame()


        counts_df_Genome = counts_df_Genome.reset_index()
        self.rewritingFile(counts_df_Genome, "countingObjectsInGenome.tsv")

        return "countingObjectsInGenome"

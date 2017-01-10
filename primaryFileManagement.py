import pandas as pa
from collections import defaultdict

def GOLabelNumberDictionnaryCreation(fileName):
    d_GOLabelToNumber = {}

    with open ("queryResults.csv", 'r') as file:
        queryResultsFile = file.read()
        queryResultsFile = queryResultsFile.replace(" ,\n", "\n")
        queryResultsFile = queryResultsFile.replace(" , ", "\t")
        queryResultsModified = open("queryResults.tsv", "w")
        queryResultsModified.write(queryResultsFile)
        queryResultsModified.close()

    queryResultsTable = pa.read_csv("queryResults.tsv", sep = "\t")

    queryResultsTable.columns = [["subject", "label", "NarrowSynonym", "BroadSynonym", "RelatedSynonym"]]

    quoteDeletion = lambda x: str(x).replace('"', '')
    queryResultsTable["subject"] = queryResultsTable["subject"].apply(quoteDeletion)

    goIsolation = lambda x: x[32:]
    queryResultsTable["subject"] = queryResultsTable["subject"].apply(goIsolation)

    spaceDeletion = lambda x: str(x).replace(" ", "")
    queryResultsTable["label"] = queryResultsTable["label"].apply(spaceDeletion)

    spaceDeletion = lambda x: str(x).replace(" ", "")
    queryResultsTable["NarrowSynonym"] = queryResultsTable["NarrowSynonym"].apply(spaceDeletion)

    spaceDeletion = lambda x: str(x).replace(" ", "")
    queryResultsTable["BroadSynonym"] = queryResultsTable["BroadSynonym"].apply(spaceDeletion)

    spaceDeletion = lambda x: str(x).replace(" ", "")
    queryResultsTable["RelatedSynonym"] = queryResultsTable["RelatedSynonym"].apply(spaceDeletion)

    d_GOLabelToNumber = dict(zip(queryResultsTable['label'], queryResultsTable['subject']))


    d_GOLabelWithSynonym = {}

    queryResultsTable = queryResultsTable.set_index(queryResultsTable['subject'])

    for index, row in queryResultsTable.iterrows():
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

    d_GOLabelWithSynonym

    return d_GOLabelToNumber, d_GOLabelWithSynonym

def translateGOTerm(listOfGOLabel, d_GOLabelToNumber):
    listOfGONumber = []

    for GOLabel in listOfGOLabel:
        if GOLabel in d_GOLabelToNumber:
            GONumber = d_GOLabelToNumber[GOLabel]
            listOfGONumber.append(GONumber)
        else:
            listOfGONumber.append(GOLabel)

    return listOfGONumber

def fixObsoleteGOTerm(listOfGOLabelAndNumber, d_GOLabelToNumber):
    listOfGONumber = []

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
            else:
                listOfGONumber.append(GOLabelAndNumber)

    return listOfGONumber

def fixWrongNorLTerm(listOfGOLabelAndNumber, d_GOLabelToNumber):
    listOfGONumber = []

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
                else:
                    listOfGONumber.append(GOLabelAndNumber)
            elif "L" in GOLabelAndNumber :
                GOLabelLWrong = GOLabelAndNumber.replace("L", "N")
                if GOLabelLWrong in d_GOLabelToNumber:
                    GONumber = d_GOLabelToNumber[GOLabelLWrong]
                    listOfGONumber.append(GONumber)
                else:
                    listOfGONumber.append(GOLabelAndNumber)
            else:
                listOfGONumber.append(GOLabelAndNumber)

    return listOfGONumber

def fixProblemsWithSynonym(listOfGOLabelAndNumber, d_GOLabelToNumber, d_GOLabelWithSynonym):
    listOfGONumber = []

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

def fixTermsIssue(listOfGOLabelAndNumber, d_GOLabelToNumber):
    listOfGONumber = []

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
            elif "al" in GOLabelAndNumber :
                GOLabelFixed = GOLabelAndNumber.replace("al", "")
                if GOLabelFixed in d_GOLabelToNumber:
                    GONumber = d_GOLabelToNumber[GOLabelFixed]
                    listOfGONumber.append(GONumber)
            elif "dependent" in GOLabelAndNumber :
                GOLabelFixed = GOLabelAndNumber.replace("dependent", "templated")
                if GOLabelFixed in d_GOLabelToNumber:
                    GONumber = d_GOLabelToNumber[GOLabelFixed]
                    listOfGONumber.append(GONumber)
            elif "sequence-specific" in GOLabelAndNumber :
                GOLabelFixed = GOLabelAndNumber.replace("sequence-specific", "") + ",sequence-specific"
                if GOLabelFixed in d_GOLabelToNumber:
                    GONumber = d_GOLabelToNumber[GOLabelFixed]
                    listOfGONumber.append(GONumber)
            else:
                listOfGONumber.append(GOLabelAndNumber)

    return listOfGONumber

def columnGOCleaning():
    resultsTable = pa.read_csv("Fichier3_Comparatif2ConditionsAvecAnnotations.csv", sep = ";")
    resultsTable = resultsTable[['Nom_Gene', 'NAME_BROAD', 'NAME_Blast2Go', 'GOs', 'EnzymeCodes', 'InterProScan']]

    resultsTable['GOs'] = resultsTable['GOs'].str.replace("C:", "")
    resultsTable['GOs'] = resultsTable['GOs'].str.replace("P:", "")
    resultsTable['GOs'] = resultsTable['GOs'].str.replace("F:", "")
    resultsTable['GOs'] = resultsTable['GOs'].str.split(";")

    d_GOLabelToNumber, d_GOLabelWithSynonym = GOLabelNumberDictionnaryCreation("queryResults.csv")

    translation = lambda x: translateGOTerm(x, d_GOLabelToNumber)
    resultsTable['GOs'] = resultsTable['GOs'].apply(translation)

    correctionProblesmSynonym = lambda x : fixProblemsWithSynonym(x,  d_GOLabelToNumber, d_GOLabelWithSynonym)
    resultsTable['GOs'] = resultsTable['GOs'].apply(correctionProblesmSynonym)

    correctionObsolete = lambda x: fixObsoleteGOTerm(x, d_GOLabelToNumber)
    resultsTable['GOs'] = resultsTable['GOs'].apply(correctionObsolete)

    correctionNorL = lambda x: fixWrongNorLTerm(x, d_GOLabelToNumber)
    resultsTable['GOs'] = resultsTable['GOs'].apply(correctionNorL)

    correctionTermsIssue = lambda x : fixTermsIssue(x, d_GOLabelToNumber)
    resultsTable['GOs'] = resultsTable['GOs'].apply(correctionTermsIssue)

    return resultsTable

def rewritingFile(newtable, fileName):
    import csv
    newtable.to_csv(fileName, "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

def main():
    resultsTable = columnGOCleaning()
    rewritingFile(resultsTable, "queryResultsGOTranslatedAndFixed.tsv")

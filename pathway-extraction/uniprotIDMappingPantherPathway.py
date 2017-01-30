import csv
import pandas as pa
import urllib2

temporaryDirectory = '../temporaryFiles/'

def httpRequestGeneOntology(url, fileName):
    req = urllib2.Request(url)
    response = urllib2.urlopen(req)
    page = response.read()
    page = page.split("\n")

    csvfile = open(temporaryDirectory + fileName + ".tsv", "w")
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow(['pathwayAccession', 'pathwayName', 'uniProtID', 'pantherSubfamilyID', 'pantherSubfamilyName'])

    for i in page[:-1]:
        columnSeparation = i.split("\t")
        writer.writerow([columnSeparation[0], columnSeparation[1], columnSeparation[4], columnSeparation[9], columnSeparation[10].replace('\r', '')])

    csvfile.close()

    return page

def cleaningFile(fileName):
    df = pa.read_csv(temporaryDirectory + fileName + ".tsv", sep = "\t")
    df['metacycPathway'] = df['metacycPathway'].str.replace("MetaCyc:", "")
    df['goLabel'] = df['goLabel'].str.replace("GO:", "")
    df.to_csv(temporaryDirectory + "output.tsv", sep= "\t", index = False, header = True, quoting = csv.QUOTE_NONE)

def main():
    fileName = 'uniprotPantherPathwayTranslation'
    test = httpRequestGeneOntology('ftp://ftp.pantherdb.org/pathway/current_release/SequenceAssociationPathway3.4.1.txt', fileName)
    cleaningFile(fileName)

main()

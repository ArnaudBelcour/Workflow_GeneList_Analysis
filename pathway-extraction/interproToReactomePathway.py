import gzip
import urllib2
from bs4 import BeautifulSoup

temporaryDirectory = '../temporaryFiles/'

def getInterproXMLCompressed(url, fileName):
    fileGZ = urllib2.urlopen(url)
    with open(temporaryDirectory + fileName + '.xml.gz','w') as interproFileCompressed:
      interproFileCompressed.write(fileGZ.read())

def uncompressInterproXML(fileCompressedName, outputFileName):
    with gzip.open(temporaryDirectory + fileName + '.xml.gz', 'r') as interproFileCompressed:
        file_content = interproFileCompressed.read()
        #interproFile = open(temporaryDirectory + fileName+ '.xml','w')
        #interproFile.write(file_content)
        #interproFile.close()
        soup = BeautifulSoup(file_content, 'lxml')
        print(soup.prettify())

def readXML(fileName):
    soup = BeautifulSoup(open(temporaryDirectory + 'interpro.xml', 'r'), 'lxml')
    print(soup.prettify())
    #with open(temporaryDirectory + 'interpro.xml','r') as interproFile:
        #print(interproFile.read())
        #soup = BeautifulSoup(interproFile.read())
        #print(soup)

def main():
    #getInterproXMLCompressed('ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', 'interpro')
    uncompressInterproXML('interpro', 'interpro')
    readXML('interpro')

#!/usr/bin/env python3

import gzip
import urllib2
from bs4 import BeautifulSoup

temporaryDirectoryDatabase = '../temporaryFiles/databases/'

def getInterproXMLCompressed(url, fileName):
    reponse = urllib2.urlopen(url, headers = {"Accept-Encoding": "gzip"})
    with GzipFile(fileobj = response) as xmlFile:
        soup = BeautifulSoup(xmlFile, 'lxml')
        print(soup.prettify())

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

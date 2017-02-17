#!/usr/bin/env python3

import gzip
import urllib2
from bs4 import BeautifulSoup

temporary_directory_database = '../temporaryFiles/databases/'

def get_interpro_xml_compressed(url, file_name):
    '''
        Requests Interpro to retrieve the interpro.xml.gz file.
        The file is analyzed to extract Interpro ID association with pathway.
    '''
    reponse = urllib2.urlopen(url, headers = {"Accept-Encoding": "gzip"})
    with GzipFile(fileobj = response) as xmlFile:
        soup = BeautifulSoup(xmlFile, 'lxml')
        print(soup.prettify())

def uncompress_interpro_xml(file_compressed_name, outputFileName):
    with gzip.open(temporaryDirectory + file_compressed_name + '.xml.gz', 'r') as interpro_file_compressed:
        file_content = interpro_file_compressed.read()
        #interproFile = open(temporaryDirectory + file_name+ '.xml','w')
        #interproFile.write(file_content)
        #interproFile.close()
        soup = BeautifulSoup(file_content, 'lxml')
        print(soup.prettify())

def read_xml(file_name):
    soup = BeautifulSoup(open(temporaryDirectory + 'interpro.xml', 'r'), 'lxml')
    print(soup.prettify())
    #with open(temporaryDirectory + 'interpro.xml','r') as interproFile:
        #print(interproFile.read())
        #soup = BeautifulSoup(interproFile.read())
        #print(soup)

def main():
    #get_interpro_xml_compressed('ftp://ftp.ebi.ac.uk/pub/databases/interpro/interpro.xml.gz', 'interpro')
    uncompress_interpro_xml('interpro', 'interpro')
    read_xml('interpro')

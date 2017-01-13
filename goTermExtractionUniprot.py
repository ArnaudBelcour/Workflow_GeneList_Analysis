import pandas as pa
import math
import requests
from bs4 import BeautifulSoup
from SPARQLWrapper import SPARQLWrapper, JSON, GET

def sparqlQueryCreation():
	df = pa.read_csv("BrowserPlasmoUniprot.tsv", skiprows = 1, sep = "\t", header=None)
	df = df[[0, 1, 2, 3]]
	df.columns = [["GenePlasmo", "PosDepart", "PosFin", "Uniprot"]]

	sparqlQueries = []

	number = 299
	listOfNumber = [0, 299]
	for iterationOf708 in range(math.ceil(len(df["Uniprot"].tolist())/300)):
		number = number + 300
		listOfNumber.append(number)
	del listOfNumber[-1]

	stringValuesProteinName = ""

	for mnemonicNameOfProtein in df["Uniprot"].tolist():
		stringValuesProteinName = stringValuesProteinName + "'" + mnemonicNameOfProtein + "' "

	for numberLimitationQuery in listOfNumber:
		stringValuesProteinName = ""

		if numberLimitationQuery != df["Uniprot"].tolist()[-1]:
			currentList = df["Uniprot"].tolist()[numberLimitationQuery:numberLimitationQuery+354]
		for mnemonicNameOfProtein in currentList:
			stringValuesProteinName = stringValuesProteinName + "'" + mnemonicNameOfProtein + "' "
		sparqlQuery = """
		PREFIX up:<http://purl.uniprot.org/core/>

		SELECT ?protein ?go
		WHERE
		{
				?protein a up:Protein .
				?protein up:mnemonic ?proteinNamed .
				?protein up:classifiedWith ?go .
				FILTER (regex(str(?go), "GO")) .
				VALUES ?proteinNamed {""" + stringValuesProteinName + """}}
		"""
		if "?proteinNamed {}" in sparqlQuery:
			continue
		else :
			sparqlQueries.append(sparqlQuery)

	return sparqlQueries

def sparqlQueryOnUniprot(sparqlQueries):
	stringValuesProteinName = 'IPO4_MOUSE'

	sparqlQuery = """
	PREFIX up:<http://purl.uniprot.org/core/>

	SELECT ?protein ?go
	WHERE
	{
			?protein a up:Protein .
			?protein up:mnemonic ?proteinNamed .
			?protein up:classifiedWith ?go .
			FILTER (regex(str(?go), "GO")) .
			VALUES ?proteinNamed {""" + stringValuesProteinName + """}}
	"""

	sparql = SPARQLWrapper("http://beta.sparql.uniprot.org/sparql")

	sparql.setQuery(sparqlQuery)
	sparql.setMethod(GET)
	sparql.setReturnFormat(JSON)

	results = sparql.query().convert()

def httprequestGOTermForAProtein(uniprotID):
    r = requests.get('http://www.uniprot.org/uniprot/' + uniprotID + '.xml')
    soup = BeautifulSoup(r.text)
    l = soup.findAll({"type", "GO"} )

    td_tag_list = soup.findAll(
                    lambda tag:tag.name == "dbreference" and
                    tag["type"] == "GO")

    GONumberList = []

    for index in range(len(td_tag_list)):
        GONumberList.append(td_tag_list[index]["id"])

    return GONumberList

def main():
	df = pa.read_csv("BrowserPlasmoUniprot.tsv", skiprows = 1, sep = "\t", header=None)
	df = df[[0, 1, 2, 3]]
	df.columns = [["GenePlasmo", "PosDepart", "PosFin", "Uniprot"]]

	GOTerms = []

	for uniprotID in df["Uniprot"].tolist():
		GOTerms.extend(requestGOTermForAProtein(uniprotID))
	
	df = pa.DataFrame(GOTerms)
	df.columns = [["GOs"]]
	df.to_csv("GOTermsPlasmoGenome.tsv", sep="\t", index = False, header = True, quoting = csv.QUOTE_NONE)
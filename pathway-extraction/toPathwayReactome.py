import requests

def httpRequestReactome(code):
    try:
        r = requests.get('http://www.reactome.org/ContentService/search/query?query=' + code +'&cluster=true')
        results = r.json()
        resultId = results['results'][0]['entries'][0]['stId']

        r = requests.get('http://www.reactome.org/ContentService/data/pathways/low/entity/' + resultId)
        print(r.json())
    except:
        print(r.json()['messages'][0])

    return r.text

def main():
    test = httpRequestReactome('0007264')

main()

import numpy as np
import pandas as pa
import unittest

import pathway_extraction.uniprot_retrieval_data as uniprot_retrieval_data

test_data_directory_uniprot = 'test_data/' + 'test_uniprot_retrieval/'

class uniprot_retrieval_data_test(unittest.TestCase):

    def test_extract_information_from_uniprot(self):
        print("\nTesting uniprot retrieval data using blast result ")
        df_data = pa.read_csv(test_data_directory_uniprot + 'data.tsv', sep='\t')
        df_data.replace(np.nan, '', regex=True, inplace=True)
        df_result = uniprot_retrieval_data.extract_information_from_uniprot(df_data)

        df_result_truth = pa.read_csv(test_data_directory_uniprot + 'result.tsv', sep='\t')

        np.testing.assert_array_equal(df_result['GOs'].tolist().sort(), df_result_truth['GOs'].tolist().sort())
        np.testing.assert_array_equal(df_result['InterProScan'].tolist(), df_result_truth['InterProScan'].tolist())

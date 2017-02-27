import numpy as np
import pandas as pa
import scipy.stats as stats
import sys
import unittest

sys.path.append("..")

from enrichmentAnalysis import EnrichmentAnalysis

test_data_directory = '../test_data/'

class enrichmentAnalysis_test(unittest.TestCase):

    def test_compute_hypergeometric_cdf(self):
        '''
        Datas are from : https://www.geneprof.org/GeneProf/tools/hypergeometric.jsp
        '''
        print("\nTesting hypergeometric cdf ")
        enrichment_analysis_test = EnrichmentAnalysis('Genes', 'data_interest_test_cdf_hypergeometric', 'data_reference_test_cdf_hypergeometric', 10, 100, 0.05, 10000)

        counts_df = pa.read_csv(test_data_directory + 'data_interest_test_cdf_hypergeometric' + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(test_data_directory + 'data_reference_test_cdf_hypergeometric' + ".tsv", sep = "\t")

        counts_df = counts_df.set_index('Genes')
        counts_df_reference = counts_df_reference.set_index('Genes')

        df_joined = counts_df.join(counts_df_reference)

        datas = [{'Genes': 'Gene_1', 'Counts': 2, 'CountsReference': 20, 'pvalue_hypergeometric': 0.6812}]
        df_joined_wih_results = pa.DataFrame(datas)

        df_joined_wih_results = df_joined_wih_results.set_index('Genes')

        for analyzed_object, row in df_joined.iterrows():
            df_joined = enrichment_analysis_test.compute_hypergeometric_test(analyzed_object, row['Counts'], row['CountsReference'], df_joined, 'under')

        np.testing.assert_array_almost_equal(df_joined['pvalue_hypergeometric'].tolist(), df_joined_wih_results['pvalue_hypergeometric'].tolist(), decimal = 4)

    def test_compute_hypergeometric_sf(self):
        '''
        Datas are from : https://fr.mathworks.com/help/stats/hygecdf.html?s_tid=gn_loc_drop
        '''
        print("\nTesting hypergeometric sf ")
        enrichment_analysis_test = EnrichmentAnalysis('Genes', 'data_interest_test_sf_hypergeometric', 'data_reference_test_sf_hypergeometric', 300, 10000, 0.05, 10000)

        counts_df = pa.read_csv(test_data_directory + 'data_interest_test_sf_hypergeometric' + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(test_data_directory + 'data_reference_test_sf_hypergeometric' + ".tsv", sep = "\t")

        counts_df = counts_df.set_index('Genes')
        counts_df_reference = counts_df_reference.set_index('Genes')

        df_joined = counts_df.join(counts_df_reference)

        datas = [{'Genes': 'Gene_1', 'Counts': 60, 'CountsReference': 2000, 'pvalue_hypergeometric': 0.5237255041}]
        df_joined_wih_results = pa.DataFrame(datas)

        df_joined_wih_results = df_joined_wih_results.set_index('Genes')

        for analyzed_object, row in df_joined.iterrows():
            df_joined = enrichment_analysis_test.compute_hypergeometric_test(analyzed_object, row['Counts'], row['CountsReference'], df_joined, 'over')

        np.testing.assert_array_almost_equal(df_joined['pvalue_hypergeometric'].tolist(), df_joined_wih_results['pvalue_hypergeometric'].tolist(), decimal = 4)

if __name__ == '__main__':
    unittest.main()

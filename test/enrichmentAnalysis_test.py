import numpy as np
import pandas as pa
import scipy.stats as stats
import sys
import unittest

sys.path.append("..")

from enrichmentAnalysis import EnrichmentAnalysis

test_data_directory = '../test_data/'


class enrichmentAnalysis_test(unittest.TestCase):

    def setUp(self):
        self.obj = EnrichmentAnalysis('Genes', 'counting_objects_in_interest', 'counting_objects_in_reference', 122, 1293, 0.05, 10000)

    def tearDown(self):
        del self.obj

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

    def test_correction_bonferroni(self):
        '''
        Datas are from : http://www.pmean.com/05/MultipleComparisons.asp
        '''
        print("\nTesting Bonferroni multiple testing correction ")
        pvalue_df = pa.read_csv(test_data_directory + 'multiple_test_data_pmean' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"

        pvalue_df = self.obj.correction_bonferroni(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory + 'multiple_test_result_pmean' + ".tsv", sep = "\t")

        np.testing.assert_array_almost_equal(pvalue_df['pValueBonferroni'].tolist(), pvalue_truth_df['PvalueBonferroni'].tolist(), decimal = 4)

    def test_correction_holm(self):
        '''
        Datas are from : http://www.pmean.com/05/MultipleComparisons.asp
        '''
        print("\nTesting Holm multiple testing correction ")
        pvalue_df = pa.read_csv(test_data_directory + 'multiple_test_data_pmean' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"

        pvalue_df = self.obj.correction_holm(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory + 'multiple_test_result_pmean' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")

        np.testing.assert_array_almost_equal(pvalue_df['pValueHolm'].tolist(), pvalue_truth_df['PvalueHolm'].tolist(), decimal = 4)

    def test_correction_benjamini_hochberg(self):
        '''
        Datas are from : http://www.pmean.com/05/MultipleComparisons.asp
        '''
        print("\nTesting Benjamini and Hochberg multiple testing correction ")
        pvalue_df = pa.read_csv(test_data_directory + 'multiple_test_data_pmean' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"

        pvalue_df = self.obj.correction_benjamini_hochberg(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory + 'multiple_test_result_pmean' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")
        np.testing.assert_array_almost_equal(pvalue_df['pValueBenjaminiHochberg'].tolist(), pvalue_truth_df['pValueBenjaminiHochberg'].tolist(), decimal = 4)

if __name__ == '__main__':
    unittest.main()

import numpy as np
import pandas as pa
import scipy.stats as stats
import unittest

from enrichmentAnalysis import EnrichmentAnalysis

test_data_directory = 'test_data/'
test_data_directory_cdf = test_data_directory + 'test_cdf/'
test_data_directory_sf = test_data_directory + 'test_sf/'
test_data_directory_multiple = test_data_directory + 'test_multiple/'

class enrichmentAnalysis_test(unittest.TestCase):

    def setUp(self):
        self.obj = EnrichmentAnalysis('Genes', 'counting_objects_in_interest', 'counting_objects_in_reference', 122, 1293, 0.05, 10000)

    def tearDown(self):
        del self.obj

    def test_percentage_calculator(self):
        print("\nTesting percentage")
        percentage_computed =  self.obj.percentage_calculator(53, 260)

        np.testing.assert_almost_equal(percentage_computed, 20.38461538461538, decimal = 5)

    def test_hypergeometric_test_on_dataframe_over(self):
        '''
        Datas are from : https://fr.mathworks.com/help/stats/hygecdf.html?s_tid=gn_loc_drop
        '''
        print("\nTesting hypergeometric sf on dataframe ")
        enrichment_analysis_test = EnrichmentAnalysis('Genes', 'data_interest_test_sf_hypergeometric', 'data_reference_test_sf_hypergeometric', 300, 10000, 0.05, 10000)

        counts_df = pa.read_csv(test_data_directory_sf + 'data_interest_test_sf_hypergeometric' + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(test_data_directory_sf + 'data_reference_test_sf_hypergeometric' + ".tsv", sep = "\t")

        counts_df.set_index('Genes', inplace = True)
        counts_df_reference.set_index('Genes', inplace = True)

        df_joined = counts_df.join(counts_df_reference)

        df_joined_wih_results = pa.read_csv(test_data_directory_sf + 'result_test_sf_hypergeometric' + ".tsv", sep = "\t")

        df_joined_wih_results.set_index('Genes', inplace = True)

        df_joined = enrichment_analysis_test.hypergeometric_test_on_dataframe(df_joined, 'over', 'CountsReference')

        np.testing.assert_array_almost_equal(df_joined['pvalue_hypergeometric'].tolist(), df_joined_wih_results['pvalue_hypergeometric'].tolist(), decimal = 4)

    def test_hypergeometric_test_on_dataframe_under(self):
        '''
        Datas are from : https://www.geneprof.org/GeneProf/tools/hypergeometric.jsp
        '''
        print("\nTesting hypergeometric cdf on dataframe ")
        enrichment_analysis_test = EnrichmentAnalysis('Genes', 'data_interest_test_cdf_hypergeometric', 'data_reference_test_cdf_hypergeometric', 10, 100, 0.05, 10000)

        counts_df = pa.read_csv(test_data_directory_cdf + 'data_interest_test_cdf_hypergeometric' + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(test_data_directory_cdf + 'data_reference_test_cdf_hypergeometric' + ".tsv", sep = "\t")

        counts_df.set_index('Genes', inplace = True)
        counts_df_reference.set_index('Genes', inplace = True)

        df_joined = counts_df.join(counts_df_reference)

        df_joined_wih_results = pa.read_csv(test_data_directory_cdf + 'result_test_cdf_hypergeometric' + ".tsv", sep = "\t")

        df_joined_wih_results.set_index('Genes', inplace = True)

        df_joined = enrichment_analysis_test.hypergeometric_test_on_dataframe(df_joined, 'under', 'CountsReference')

        np.testing.assert_array_almost_equal(df_joined['pvalue_hypergeometric'].tolist(), df_joined_wih_results['pvalue_hypergeometric'].tolist(), decimal = 4)

    def test_correction_bonferroni(self):
        '''
        Datas are from : http://www.pmean.com/05/MultipleComparisons.asp
        '''
        print("\nTesting Bonferroni multiple testing correction ")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_pmean' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"

        pvalue_df = self.obj.correction_bonferroni(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_pmean' + ".tsv", sep = "\t")

        np.testing.assert_array_almost_equal(pvalue_df['pValueBonferroni'].tolist(), pvalue_truth_df['PvalueBonferroni'].tolist(), decimal = 4)

    def test_correction_holm(self):
        '''
        Datas are from : http://www.pmean.com/05/MultipleComparisons.asp
        '''
        print("\nTesting Holm multiple testing correction ")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_pmean' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"

        pvalue_df = self.obj.correction_holm(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_pmean' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")

        np.testing.assert_array_almost_equal(pvalue_df['pValueHolm'].tolist(), pvalue_truth_df['PvalueHolm'].tolist(), decimal = 4)

    def test_correction_benjamini_hochberg(self):
        '''
        Datas are from : http://www.pmean.com/05/MultipleComparisons.asp
                        https://journal.r-project.org/archive/2014-2/conde-alvarez.pdf
        '''
        print("\nTesting Benjamini and Hochberg multiple testing correction ")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_pmean' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"

        pvalue_df = self.obj.correction_benjamini_hochberg(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_pmean' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")
        np.testing.assert_array_almost_equal(pvalue_df['pValueBenjaminiHochberg'].tolist(), pvalue_truth_df['pValueBenjaminiHochberg'].tolist(), decimal = 4)

    def test_correction_sgof_G(self):
        '''
        Datas are from : http://acraaj.webs.uvigo.es/SGoFReadme.htm
        '''
        print("\nTesting SGoF multiple testing correction using G test")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_sgof_G_test' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"
        self.obj.object_to_analyze= "pvalue_hypergeometric"
        pvalue_df = self.obj.correction_sgof(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_sgof_G_test' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")

        np.testing.assert_array_equal(pvalue_df['pValueSGoF'].tolist(), pvalue_truth_df['pValueSGoFG'].tolist())

    def test_correction_sgof_bino(self):
        '''
        Datas are from : http://acraaj.webs.uvigo.es/SGoFReadme.htm
        '''
        print("\nTesting SGoF multiple testing correction using binomial test ")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_sgof_binomial' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"
        self.obj.object_to_analyze= "pvalue_hypergeometric"
        pvalue_df = self.obj.correction_sgof(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_sgof_binomial' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")

        np.testing.assert_array_equal(pvalue_df['pValueSGoF'].tolist(), pvalue_truth_df['pValueSGoFBino'].tolist())

if __name__ == '__main__':
    unittest.main()

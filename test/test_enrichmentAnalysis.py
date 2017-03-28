import numpy as np
import os
import pandas as pa
import scipy.stats as stats
import unittest

from fileManagement import FileManagement
from enrichmentAnalysis import EnrichmentAnalysis, GOEnrichmentAnalysis
from unittest.mock import patch

test_data_directory = 'test_data/'
test_data_directory_cdf = test_data_directory + 'test_cdf/'
test_data_directory_sf = test_data_directory + 'test_sf/'
test_data_directory_multiple = test_data_directory + 'test_multiple/'
test_data_directory_enrichment = test_data_directory + 'test_enrichment/'

class enrichmentAnalysis_test(unittest.TestCase):

    def setUp(self):
        self.obj = EnrichmentAnalysis('Genes', 'counting_objects_in_interest', 'counting_objects_in_reference', 122, 1293, 0.05, 10000)

    def tearDown(self):
        del self.obj

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

    def test_normal_approximation_on_dataframe_over(self):
        '''
        Datas have been invented for the test.
        Results for the hypergeometric test are from : https://www.geneprof.org/GeneProf/tools/hypergeometric.jsp
        pvalue_hypergeometric found : 4.834533775884863e-8
        '''
        print("\nTesting normal approximation sf on dataframe ")
        enrichment_analysis_test = EnrichmentAnalysis('Genes', 'data_interest_test_sf_hypergeometric', 'data_reference_test_sf_hypergeometric', 10000, 100000, 0.05, 10000)

        counts_df = pa.read_csv(test_data_directory_sf + 'data_interest_test_sf_hypergeometric_normal' + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(test_data_directory_sf + 'data_reference_test_sf_hypergeometric_normal' + ".tsv", sep = "\t")

        counts_df.set_index('Genes', inplace = True)
        counts_df_reference.set_index('Genes', inplace = True)

        df_joined = counts_df.join(counts_df_reference)

        df_joined['pvalue_normal_approximation'] = df_joined.apply(enrichment_analysis_test.compute_normal_approximation, args=('CountsReference', 'over'), axis=1)
        df_joined = enrichment_analysis_test.hypergeometric_test_on_dataframe(df_joined, 'over', 'CountsReference')

        np.testing.assert_array_almost_equal(df_joined['pvalue_hypergeometric'].tolist(), df_joined['pvalue_normal_approximation'].tolist(), decimal = 4)

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
        Datas and results are from : www.biostathandbook.com/multiplecomparisons.html
        Data are from the article : onlinelibrary.wiley.com/doi/10.1002/ijc.28513/full
        '''
        print("\nTesting Benjamini and Hochberg multiple testing correction ")

        df_data = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_BH.tsv', sep='\t')
        self.obj.statistic_method = "pvalue_hypergeometric"

        df_data = self.obj.correction_benjamini_hochberg(df_data)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_BH.tsv', sep='\t')

        np.testing.assert_array_almost_equal(df_data['pValueBenjaminiHochberg'].tolist(), pvalue_truth_df['pValueBenjaminiHochberg'].tolist(), decimal = 4)

    def test_correction_sgof_G(self):
        '''
        Datas have been created for the example.
        To obtain the results, they have been used on MATLAB script for SGoF here : http://acraaj.webs.uvigo.es/software/matlab_sgof.m
        '''
        print("\nTesting SGoF multiple testing correction using G test")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_sgof_G_test' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"
        self.obj.object_to_analyze= "pvalue_hypergeometric"
        pvalue_df = self.obj.correction_sgof(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_sgof_G_test' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")

        np.testing.assert_array_equal(pvalue_df['pValueSGoF'].tolist(), pvalue_truth_df['pValueSGoF'].tolist())

    def test_correction_sgof_bino(self):
        '''
        Datas have been created for the example.
        To obtain the results, they have been used on MATLAB script for SGoF here : http://acraaj.webs.uvigo.es/software/matlab_sgof.m
        '''
        print("\nTesting SGoF multiple testing correction using binomial test ")
        pvalue_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_data_sgof_binomial' + ".tsv", sep = "\t")

        self.obj.statistic_method = "pvalue_hypergeometric"
        self.obj.object_to_analyze= "pvalue_hypergeometric"
        pvalue_df = self.obj.correction_sgof(pvalue_df)

        pvalue_truth_df = pa.read_csv(test_data_directory_multiple + 'multiple_test_result_sgof_binomial' + ".tsv", sep = "\t")
        pvalue_truth_df = pvalue_truth_df.sort_values(by = "pvalue_hypergeometric")

        np.testing.assert_array_equal(pvalue_df['pValueSGoF'].tolist(), pvalue_truth_df['pValueSGoF'].tolist())

    def test_error_rate_adjustement_bonferroni(self):
        '''
        Datas and results are from : www.biostathandbook.com/multiplecomparisons.html
        '''
        print("\nTesting error rate adjustement Bonferroni ")
        datas = {'pvalue_hypergeometric':[0.001,0.008,0.039,0.041,0.042,0.06,0.074,0.205,0.212,0.216,0.222,
                                    0.251,0.269,0.275,0.34,0.341,0.384,0.569,0.594,0.696,0.762,0.94,0.942,0.975,0.986]}
        df = pa.DataFrame(datas)
        error_rate_adjusted = self.obj.error_rate_adjustement_bonferroni(df)

        self.assertEqual(error_rate_adjusted, 0.002)

    def test_error_rate_adjustement_sidak(self):
        '''
        Datas have been created for the example (the only important thing here is the numver of pvalue).
        The example and the result are from : www.spc.univ-lyon1.fr/polycop/comparaisons multiples.htm
        '''
        print("\nTesting error rate adjustement Sidak ")
        datas_10 = {'pvalue_hypergeometric_10':[0.01,0.02,0.3,0.02,0.05,0.07,0.9,0.001,0.09,0.008]}
        df_10_pvalue = pa.DataFrame(datas_10)
        error_rate_adjusted_10 = self.obj.error_rate_adjustement_sidak(df_10_pvalue)

        datas_20 = {'pvalue_hypergeometric_20':[0.01,0.02,0.05,0.04,0.2,0.04,0.9,0.05,0.06,0.0545,
                                                0.048766,0.02,0.04,0.03,0.365,0.21,0.0234,0.2,0.156]}
        df_20_pvalue = pa.DataFrame(datas_20)
        error_rate_adjusted_20 = self.obj.error_rate_adjustement_sidak(df_20_pvalue)

        np.testing.assert_array_almost_equal([error_rate_adjusted_10, error_rate_adjusted_20], [0.0051, 0.0026], decimal = 4)

    def test_enrichment_analysis(self):
        '''
        Datas are from an enrichment analysis course.
        Use a mock to simulate yes_or_no input for approximation and change global variable.
        '''
        print("\nTesting enrichment analysis ")
        file_management = FileManagement(test_data_directory_enrichment + 'counting_objects_in_genome.tsv')

        with patch('fileManagement.temporary_directory', test_data_directory_enrichment):
            go_number_go_labels = file_management.go_label_number_dictionary_creation_from_http(specification='inverse')

        with patch('builtins.input', return_value='n'):
            with patch('enrichmentAnalysis.temporary_directory', test_data_directory_enrichment):
                with patch('enrichmentAnalysis.output_directory', test_data_directory_enrichment):
                    go_enrichment_analysis = GOEnrichmentAnalysis('GOs', 'counting_objects_in_interest', 'counting_objects_in_genome',
                                                                    11, 38660, 0.05, 10000, go_number_go_labels)
                    go_enrichment_analysis.enrichment_analysis()

                    results = pa.read_csv(test_data_directory_enrichment + 'pValuesOfGOs_over.tsv', sep='\t', float_precision='high', skiprows=1)
                    results_truth = pa.read_csv(test_data_directory_enrichment + 'overRepresentation_genesSet1.tsv', sep='\t')
                    results.sort_values('GOs', inplace=True)
                    results_truth.sort_values('GOs', inplace=True)

                    np.testing.assert_array_almost_equal(results['pvalue_hypergeometric'].tolist(), results_truth['pvalue_hypergeometric'].tolist(), decimal = 4)
                    np.testing.assert_array_equal(results['GOLabel'].tolist(), results_truth['Labels'].tolist())

                    os.remove(test_data_directory_enrichment + 'pValuesOfGOs_over.tsv')
                    os.remove(test_data_directory_enrichment + 'pValuesOfGOs_under.tsv')
                    os.remove(test_data_directory_enrichment + 'go_number_label.tsv')
                    os.remove(test_data_directory_enrichment + 'significativesGOs_over.tsv')
                    os.remove(test_data_directory_enrichment + 'significativesGOs_under.tsv')

if __name__ == '__main__':
    unittest.main()

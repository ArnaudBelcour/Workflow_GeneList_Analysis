import pandas as pa
import scipy.stats as stats
import sys
import unittest

from pandas.util.testing import assert_frame_equal
sys.path.append("..")

from enrichmentAnalysis import EnrichmentAnalysis

test_data_directory = '../test_data/'

class enrichmentAnalysis_test(unittest.TestCase):

    def test_compute_hypergeometric(self):
        enrichment_analysis_test = EnrichmentAnalysis('Animals', 'counting_objects_in_interest', 'counting_objects_in_reference', 122, 1293, 0.05, 10000)

        counts_df = pa.read_csv(test_data_directory + 'counting_objects_in_interest' + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(test_data_directory + 'counting_objects_in_reference' + ".tsv", sep = "\t")
 
        counts_df = counts_df.set_index('Animals')
        counts_df_reference = counts_df_reference.set_index('Animals')

        df_joined = counts_df.join(counts_df_reference)

        df_joined_wih_results = df_joined.copy()

        for analyzed_object, row in df_joined.iterrows():
            df_joined = enrichment_analysis_test.compute_hypergeometric_test(analyzed_object, row['Counts'], row['CountsReference'], df_joined, 'over')
            df_joined_wih_results.set_value(analyzed_object, 'pvalue_hypergeometric', stats.hypergeom.sf(row['Counts'] - 1, 1293, row['CountsReference'], 122))

        comparisons = df_joined.values == df_joined_wih_results.values

        self.assertEqual(all(comparisons.tolist()), True)

if __name__ == '__main__':
    unittest.main()
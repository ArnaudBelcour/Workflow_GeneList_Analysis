#!/usr/bin/env python3

import csv
import math
import pandas as pa
import scipy.stats as stats
import six

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
output_directory = 'outputFiles/'

class EnrichmentAnalysis():

    '''
        Performs an enrichment analysis using an hypergeometric test (also known as Fisher's exact test) and multiple correction testing.
        To do this you need to enter some values (using the set_something() function) :
            -file of interest : the name of your file (with the extension) containing the occurrences of each objects from a sample you want to analyze.
            -file of reference : the name of your file (with the extentions) containing the occurrences of each objects from a population.
            -number of analyzed object of interest : the number of objects in your sample (for example the number of differentially expressed genes in a list).
            -number of analyzed object of reference : the number of objects in your population (for example the number of genes in the genome of your species).
            -alpha : the alpha threshold also known as type I error.
            -normal approximation threshold : the threshold separating the hypergeometric test (which runs very slowly when using big numbers) and normal approximation.
    '''

    def __init__(self, column_name, file_of_interest_name, file_of_reference_name, number_of_object_of_interest, number_of_genes_in_reference, alpha, threshold_normal_approximation):
        self._object_to_analyze = column_name
        self._output_columns = ['Counts', 'CountsReference', 'Percentage' + self.object_to_analyze + 'InInterest', 'Percentage' + self.object_to_analyze + 'InReference', 'pvalue_hypergeometric', 'pValueBonferroni', 'pValueHolm', 'pValueSGoF', 'pValueBenjaminiHochberg']
        self._file_of_interest = file_of_interest_name
        self._file_of_reference = file_of_reference_name
        self._number_of_analyzed_object_of_interest = number_of_object_of_interest
        self._number_of_analyzed_object_of_reference = number_of_genes_in_reference
        self._alpha = alpha
        self._normal_approximation_threshold = threshold_normal_approximation
        self._statistic_method = ""
        self.multiple_test_names = ['Sidak', 'Bonferroni', 'Holm', 'SGoF', 'BenjaminiHochberg']

    @property
    def object_to_analyze(self):
         return self._object_to_analyze

    @object_to_analyze.setter
    def object_to_analyze(self, object_name):
         self._object_to_analyze = object_name

    @property
    def output_columns(self):
        return self._output_columns

    @output_columns.setter
    def output_columns(self, index, column_name):
        self._output_columns[index] = column_name

    @property
    def file_of_interest(self):
        return self._file_of_interest

    @file_of_interest.setter
    def file_of_interest(self, fileName):
        self._file_of_interest = fileName

    @property
    def file_of_reference(self):
        return self._file_of_reference

    @file_of_reference.setter
    def file_of_reference(self, fileName):
        self._file_of_reference = fileName

    @property
    def number_of_analyzed_object_of_interest(self):
        return self._number_of_analyzed_object_of_interest

    @number_of_analyzed_object_of_interest.setter
    def number_of_analyzed_object_of_interest(self, value):
        if value > self.number_of_analyzed_object_of_reference:
            raise ValueError("The number of objects in your sample of interest is greater than the number of objects in the reference.")
        else:
            self._number_of_analyzed_object_of_interest = value

    @property
    def number_of_analyzed_object_of_reference(self):
        return self._number_of_analyzed_object_of_reference

    @number_of_analyzed_object_of_reference.setter
    def number_of_analyzed_object_of_reference(self, value):
        if value < self.number_of_analyzed_object_of_interest:
            raise ValueError("The number of objects in the reference is smaller than the number of objects in your sample of interest.")
        else:
            self._number_of_analyzed_object_of_reference = value

    @property
    def alpha(self):
        return self._alpha

    @alpha.setter
    def alpha(self, value):
        self._alpha = value

    @property
    def statistic_method(self):
        return self._statistic_method

    @statistic_method.setter
    def statistic_method(self, method_name):
        self._statistic_method = method_name

    @property
    def normal_approximation_threshold(self):
        return self._normal_approximation_threshold

    @normal_approximation_threshold.setter
    def normal_approximation_threshold(self, value):
        self._normal_approximation_threshold = value

    def hypergeometric_test_on_dataframe(self, df, over_or_underrepresentation, reference_column):
        analyzed_objects_with_hypergeo_test_nan = []

        approximation_threshold = self.normal_approximation_threshold

        for analyzed_object, row in df.iterrows():
            if math.isnan(df.get_value(analyzed_object, reference_column)):
                df = df.drop([analyzed_object])
            else:
                if row['Counts'] < approximation_threshold:
                    self.compute_hypergeometric_test(analyzed_object, row['Counts'], row[reference_column], df, over_or_underrepresentation)

                elif row['Counts'] > approximation_threshold:
                    self.compute_normal_approximation(analyzed_object, row['Counts'], row[reference_column], df, over_or_underrepresentation)

                if math.isnan(df.get_value(analyzed_object, self.statistic_method)):
                    analyzed_objects_with_hypergeo_test_nan.append(analyzed_object)
                    df = df.drop([analyzed_object])
                df = df.sort_values(self.statistic_method)

        return df

    def compute_hypergeometric_test(self, analyzed_object, number_of_object_in_interest, number_of_object_in_reference, df, over_or_underrepresentation):
        if over_or_underrepresentation == "over":
            pvalue_hypergeo = stats.hypergeom.sf(number_of_object_in_interest - 1, self.number_of_analyzed_object_of_reference, number_of_object_in_reference, self.number_of_analyzed_object_of_interest)

        if over_or_underrepresentation == "under":
            pvalue_hypergeo = stats.hypergeom.cdf(number_of_object_in_interest, self.number_of_analyzed_object_of_reference, number_of_object_in_reference, self.number_of_analyzed_object_of_interest)

        df.set_value(analyzed_object, 'pvalue_hypergeometric', pvalue_hypergeo)
        self.statistic_method = "pvalue_hypergeometric"

        return df

    def compute_normal_approximation(self, analyzedObject, number_of_object_in_interest, number_of_object_in_reference, df, over_or_underrepresentation):
        p = number_of_object_in_reference / self.number_of_analyzed_object_of_reference
        q = 1 - p
        t = self.number_of_analyzed_object_of_interest / self.number_of_analyzed_object_of_reference

        mu = self.number_of_analyzed_object_of_interest  * p
        sigma = math.sqrt(self.number_of_analyzed_object_of_interest  * p * q * (1 - t))

        if over_or_underrepresentation == "over":
            pValueNormal = stats.norm.sf(number_of_object_in_interest + 1, loc = mu, scale = sigma)

        if over_or_underrepresentation == "under":
            pValueNormal = stats.norm.cdf(number_of_object_in_interest, loc = mu, scale = sigma)

        df.set_value(analyzedObject, 'pvalue_normal_approximation', pValueNormal)
        self.output_columns[4] = 'pvalue_normal_approximation'
        self.statistic_method = 'pvalue_normal_approximation'

        return df

    def counting_approximation(self, df):
        for analyzed_object, row in df.iterrows():
            df.set_value(analyzed_object, 'CountsTotal', row['Counts'] + row['CountsReference'])

        return df

    def percentage_calculation(self, numerator, denominator):
        percentage = (numerator / denominator) * 100

        return percentage

    def multiple_testing_correction(self, df):
        df = df.sort_values([self.statistic_method])

        df = self.correction_bonferroni(df)
        df = self.correction_benjamini_hochberg(df)
        df = self.correction_holm(df)
        df = self.correction_sgof(df)

        significative_objects = {}

        for multiple_test_name in self.multiple_test_names:
            if multiple_test_name == 'Sidak':
                error_rate = self.error_rate_adjustement_sidak(df)
            elif multiple_test_name == 'Bonferroni':
                error_rate = self.error_rate_adjustement_bonferroni(df)
            if multiple_test_name in ['Sidak', 'Bonferroni']:
                object_significatives = self.selection_object_with_adjusted_error_rate(error_rate, df)
            elif multiple_test_name in ['Holm', 'BenjaminiHochberg']:
                object_significatives = self.selection_object_with_adjusted_pvalue(multiple_test_name, df)
            elif multiple_test_name == 'SGoF':
                object_significatives = self.selection_object_with_sgof(multiple_test_name, df)

            significative_objects[multiple_test_name] = object_significatives

        return df, significative_objects

    def writing_output(self, df, significative_objects, over_or_underrepresentation, approximation_yes_or_no, yes_answers):
        '''
        For the second results file (file with significative objects):
        Results are written using sorted(dictionnary), so the list of result corresponds to : Sidak (position 4 in the list), Bonferroni (position 1),
        Holm (position 2), SGoF (position 3) and Benjamini and Hochberg (position 0).
        '''
        df = df.sort_values(['pValueBenjaminiHochberg'])

        if approximation_yes_or_no in yes_answers:
            self.output_columns[1] = "CountsTotal"
            df = df[self.output_columns]
        else:
            df = df[self.output_columns]

        if over_or_underrepresentation == 'over':
            df.to_csv(output_directory + "pValuesOf" + self.object_to_analyze + "_over.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)
        elif over_or_underrepresentation == 'under':
            df.to_csv(output_directory + "pValuesOf" + self.object_to_analyze + "_under.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)

        if over_or_underrepresentation == 'over':
            csvfile = open(output_directory + "significatives" + self.object_to_analyze + "_over.tsv", "w", newline = "")
        elif over_or_underrepresentation == 'under':
            csvfile = open(output_directory + "significatives" + self.object_to_analyze + "_under.tsv", "w", newline = "")

        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow([self.object_to_analyze + 'Sidak', self.object_to_analyze + 'Bonferroni', self.object_to_analyze + 'Holm', self.object_to_analyze + 'SGoF', self.object_to_analyze + 'BenjaminiHochberg'])

        number_significatives_per_method = {}

        for method in significative_objects:
            number_significatives_per_method[method] = len(significative_objects[method])

        max_significatives_method = max(number_significatives_per_method, key = number_significatives_per_method.get)

        for index in range(len(significative_objects[max_significatives_method])):
            results = []
            for method in sorted(significative_objects):
                if index in range(len(significative_objects[method])):
                    object_significatives_value =  significative_objects[method][index]
                else :
                    object_significatives_value =  'nan'
                results.append(object_significatives_value)

            writer.writerow([results[4], results[1], results[2], \
            results[3], results[0]])

        csvfile.close()

    def correction_bonferroni(self, df):

        def pvalue_correction_bonferroni(x) :
            number_of_test = len(df.index)
            q_value = x * number_of_test
            if q_value > 1 :
                q_value = 1

            return q_value

        df['pValueBonferroni'] = df[self.statistic_method].apply(pvalue_correction_bonferroni)

        return df

    def correction_benjamini_hochberg(self, df):
        df = df.sort_values(by = self.statistic_method, ascending = True)

        number_of_test = len(df.index)

        for analyzed_object, row in df.iterrows():
            rank = df.index.get_loc(analyzed_object)
            pvalue_correction_benjamini_hochberg = row[self.statistic_method] * (number_of_test / (rank + 1))

            if pvalue_correction_benjamini_hochberg > 1:
                pvalue_correction_benjamini_hochberg = 1
            df.set_value(analyzed_object, 'pValueBenjaminiHochberg', pvalue_correction_benjamini_hochberg)

        return df

    def correction_holm(self, df):
        df = df.sort_values(by = self.statistic_method, ascending = True)

        number_of_test = len(df.index)
        pvalue_max = 0

        for analyzed_object, row in df.iterrows():
            rank = df.index.get_loc(analyzed_object)
            pvalue_correction_holm = row[self.statistic_method] * (number_of_test - rank)
            if pvalue_correction_holm > 1:
                pvalue_correction_holm = 1
            if pvalue_correction_holm > pvalue_max:
                pvalue_max = pvalue_correction_holm
            if pvalue_max > pvalue_correction_holm:
                pvalue_correction_holm = pvalue_max
            if pvalue_correction_holm > self.alpha and pvalue_max < pvalue_correction_holm:
                pvalue_max = pvalue_correction_holm

            df.set_value(analyzed_object, 'pValueHolm', pvalue_correction_holm)

        return df

    def correction_sgof(self, df):
        F = len(df.index) * self.alpha
        df = df.sort_values("pvalue_hypergeometric")
        R = 0

        R = (df['pvalue_hypergeometric'] < 0.05).sum()

        df = df.reset_index()

        object_significatives = []
        row_number = 0
        while stats.binom_test(R, len(df), p = self.alpha) < self.alpha and R > 0:
            df.set_value(row_number, 'pValueSGoF', 'significant')
            object_significatives.append(df.iloc[row_number][self.object_to_analyze])
            R = R - 1
            row_number = row_number + 1

        df = df.set_index(self.object_to_analyze)

        for analyzed_object, row in df.iterrows():
            try:
                if analyzed_object not in object_significatives:
                    df.set_value(analyzed_object, 'pValueSGoF', 'nonSignificant')
            except:
                df.set_value(analyzed_object, 'pValueSGoF', 'nonSignificant')

        return df

    def error_rate_adjustement_bonferroni(self, df):
        error_rate_adjusted = self.alpha / len(df.index)

        return error_rate_adjusted

    def error_rate_adjustement_sidak(self, df):
        error_rate_adjusted = (1 - math.pow((1 - self.alpha), (1 / len(df.index))))

        return error_rate_adjusted

    def selection_object_with_adjusted_error_rate(self, error_rate, df):
        object_significatives = []

        for analyzed_object, row in df.iterrows():
            if row[self.statistic_method] < error_rate :
                object_significatives.append(analyzed_object)

        return object_significatives

    def selection_object_with_adjusted_pvalue(self, method_name, df):
        object_significatives = []

        for analyzed_object, row in df.iterrows():
            if row['pValue' + method_name] < self.alpha:
                object_significatives.append(analyzed_object)

        return object_significatives

    def selection_object_with_sgof(self, method_name, df):
        object_significatives = []

        for analyzed_object, row in df.iterrows():
            if row['pValue' + method_name] == 'significant' :
                object_significatives.append(analyzed_object)

        return object_significatives

    def enrichment_analysis(self):

        counts_df = pa.read_csv(temporary_directory + self.file_of_interest + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(temporary_directory + self.file_of_reference + ".tsv", sep = "\t")

        counts_df = counts_df.set_index(self.object_to_analyze)
        counts_df_reference = counts_df_reference.set_index(self.object_to_analyze)

        df_joined = counts_df.join(counts_df_reference)

        yes_answers = ['yes', 'y', 'oui', 'o']
        yes_or_no = input("Is this an approximation of the reference? ")

        for analyzed_object, row in df_joined.iterrows():
            df_joined.set_value(analyzed_object, 'Percentage' + self.object_to_analyze + 'InInterest', self.percentage_calculation(row['Counts'], self.number_of_analyzed_object_of_interest))

        if yes_or_no in yes_answers:
            df_joined = self.counting_approximation(df_joined)
            count_column_name = 'CountsTotal'
        else:
            count_column_name = 'CountsReference'

        for analyzed_object, row in df_joined.iterrows():
            df_joined.set_value(analyzed_object, 'Percentage' + self.object_to_analyze + 'InReference', self.percentage_calculation(row[count_column_name], self.number_of_analyzed_object_of_reference))

        over_unders = ['over', 'under']
        for over_under in over_unders:
            df_joined_over_under = self.hypergeometric_test_on_dataframe(df_joined, over_under, count_column_name)
            df_joined_over_under, significative_objects = self.multiple_testing_correction(df_joined_over_under)
            self.writing_output(df_joined_over_under, significative_objects, over_under, yes_or_no, yes_answers)


class GOEnrichmentAnalysis(EnrichmentAnalysis):

    def __init__(self, column_name, file_of_interest_name, file_of_reference_name, number_of_object_of_interest, number_of_genes_in_reference, alpha, threshold_normal_approximation, d_go_label_to_number):
        EnrichmentAnalysis.__init__(self, column_name, file_of_interest_name, file_of_reference_name, number_of_object_of_interest, number_of_genes_in_reference, alpha, threshold_normal_approximation)
        self.output_columns.append("GOLabel")
        self._gos_labels_to_numbers = d_go_label_to_number

    @property
    def gos_labels_to_numbers(self):
        return self._gos_labels_to_numbers

    @gos_labels_to_numbers.setter
    def gos_labels_to_numbers(self, go_dictionnary):
        self._gos_labels_to_numbers = go_dictionnary

    def tranlsation_go_number_to_go_label(self, go_numbers, d_go_label_to_number):
        go_labels = []

        for go_number in go_numbers:
            if go_number in d_go_label_to_number:
                go_labels.append(d_go_label_to_number[go_number])

        return go_labels

    def multiple_testing_correction(self, df):
        df = df.sort_values([self.statistic_method])

        df = self.correction_bonferroni(df)
        df = self.correction_benjamini_hochberg(df)
        df = self.correction_holm(df)
        df = self.correction_sgof(df)

        significative_objects = {}
        translation_gos_labels_to_numbers = self.gos_labels_to_numbers

        for multiple_test_name in self.multiple_test_names:
            if multiple_test_name == 'Sidak':
                error_rate = self.error_rate_adjustement_sidak(df)
            elif multiple_test_name == 'Bonferroni':
                error_rate = self.error_rate_adjustement_bonferroni(df)
            if multiple_test_name in ['Sidak', 'Bonferroni']:
                object_significatives = self.selection_object_with_adjusted_error_rate(error_rate, df)
            elif multiple_test_name in ['Holm', 'BenjaminiHochberg']:
                object_significatives = self.selection_object_with_adjusted_pvalue(multiple_test_name, df)
            elif multiple_test_name == 'SGoF':
                object_significatives = self.selection_object_with_sgof(multiple_test_name, df)

            go_label_significatives = self.tranlsation_go_number_to_go_label(object_significatives, translation_gos_labels_to_numbers)
            significative_objects[multiple_test_name] = go_label_significatives

        for go, row in df.iterrows():
            if go in translation_gos_labels_to_numbers:
                df.set_value(go, 'GOLabel', translation_gos_labels_to_numbers[go])

        return df, significative_objects

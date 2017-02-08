#!/usr/bin/env python3

import csv
import math
import pandas as pa
import scipy.stats as stats
import six
import sys
from ast import literal_eval

input_directory = "inputFiles/"
temporary_directory = 'temporaryFiles/'
output_directory = 'outputFiles/'

class EnrichmentAnalysis():

    def __init__(self, column_name):
        self.object_to_analyze = column_name
        self.output_columns = ['Counts', 'CountsReference', 'Percentage' + self.get_object_to_analyze() + 'InInterest', 'Percentage' + self.get_object_to_analyze() + 'InReference', 'pvalue_hypergeometric', 'pValueBonferroni', 'pValueHolm', 'pValueSGoF', 'pValueBenjaminiHochberg']
        self.file_of_interest = ""
        self.file_of_reference = ""
        self.number_of_analyzed_object_of_interest = 0
        self.number_of_analyzed_object_of_reference = 0
        self.alpha = 0.00
        self.statistic_method = ""

    def get_object_to_analyze(self):
        return self.object_to_analyze

    def get_output_columns(self):
        return self.output_columns

    def get_file_of_interest(self):
        return self.file_of_interest

    def get_file_of_reference(self):
        return self.file_of_reference

    def get_number_of_analyzed_object_of_interest(self):
        return self.number_of_analyzed_object_of_interest

    def get_number_of_analyzed_object_of_reference(self):
        return self.number_of_analyzed_object_of_reference

    def get_alpha(self):
        return self.alpha

    def get_statistic_method(self):
        return self.statistic_method

    def set_file_of_interest(self, fileName):
        self.file_of_interest = fileName

    def set_file_of_reference(self, fileName):
        self.file_of_reference = fileName

    def set_output_columns(self, index, column_name):
        self.output_columns[index] = column_name

    def set_number_of_analyzed_object_of_interest(self, value):
        self.number_of_analyzed_object_of_interest = value

    def set_number_of_analyzed_object_of_reference(self, value):
        self.number_of_analyzed_object_of_reference = value

    def set_alpha(self, value):
        self.alpha = value

    def set_statistic_method(self, method_name):
        self.statistic_method = method_name

    def hypergeometric_test_on_dataframe(self, df, over_or_underrepresentation, genome_columns):
        analyzed_objects_with_hypergeo_test_nan = []

        for analyzed_object, row in df.iterrows():
            if math.isnan(df.get_value(analyzed_object, genome_columns)):
                df = df.drop([analyzed_object])
            else:
                if row['Counts'] < 10000:
                    self.compute_hypergeometric_test(analyzed_object, row['Counts'], row[genome_columns], df, over_or_underrepresentation)

                elif row['Counts'] > 10000:
                    self.compute_normal_approximation(analyzed_object, row['Counts'], row[genome_columns], df, over_or_underrepresentation)

                if math.isnan(df.get_value(analyzed_object, self.get_statistic_method())):
                    analyzed_objects_with_hypergeo_test_nan.append(analyzed_object)
                    df = df.drop([analyzed_object])
                df = df.sort_values(self.get_statistic_method())

        return df

    def compute_hypergeometric_test(self, analyzed_object, number_of_object_in_interest, number_of_object_in_reference, df, over_or_underrepresentation):
        if over_or_underrepresentation == "over":
            pvalue_hypergeo = stats.hypergeom.sf(number_of_object_in_interest - 1, self.get_number_of_analyzed_object_of_reference(), number_of_object_in_reference, self.get_number_of_analyzed_object_of_interest())

        if over_or_underrepresentation == "under":
            pvalue_hypergeo = stats.hypergeom.cdf(number_of_object_in_interest, self.get_number_of_analyzed_object_of_reference(), number_of_object_in_reference, self.get_number_of_analyzed_object_of_interest())

        df.set_value(analyzed_object, 'pvalue_hypergeometric', pvalue_hypergeo)
        self.set_statistic_method("pvalue_hypergeometric")

        return df

    def compute_normal_approximation(self, analyzedObject, number_of_object_in_interest, number_of_object_in_reference, df, over_or_underrepresentation):
        p = number_of_object_in_reference / self.get_number_of_analyzed_object_of_reference()
        q = 1 - p
        t = self.get_number_of_analyzed_object_of_interest() / self.get_number_of_analyzed_object_of_reference()

        mu = self.get_number_of_analyzed_object_of_interest()  * p
        sigma = math.sqrt(self.get_number_of_analyzed_object_of_interest()  * p * q * (1 - t))

        if over_or_underrepresentation == "over":
            pValueNormal = stats.norm.sf(number_of_object_in_interest + 1, loc = mu, scale = sigma)

        if over_or_underrepresentation == "under":
            pValueNormal = stats.norm.cdf(number_of_object_in_interest, loc = mu, scale = sigma)

        df.set_value(analyzedObject, 'pvalue_normal_approximation', pValueNormal)
        self.set_output_columns(4, 'pvalue_normal_approximation')
        self.set_statistic_method('pvalue_normal_approximation')

        return df

    def counting_approximation(self, df):
        for analyzed_object, row in df.iterrows():
            df.set_value(analyzed_object, 'CountsTotal', row['Counts'] + row['CountsReference'])

        return df

    def percentage_calculation(self, numerator, denominator):
        percentage = (numerator / denominator) * 100

        return percentage

    def multiple_testing_correction(self, df):
        df = df.sort_values([self.get_statistic_method()])

        df = self.correction_bonferroni(df)
        df = self.correction_benjamini_hochberg(df)
        df = self.correction_holm(df)
        df = self.correction_sgof(df)

        significative_objects = {}

        error_rate_sidak = self.error_rate_adjustement_sidak(df)
        object_significatives_Sidak = self.selection_object_with_adjusted_error_rate(error_rate_sidak, df)
        significative_objects['Sidak'] = object_significatives_Sidak

        error_rate_bonferroni = self.error_rate_adjustement_bonferroni(df)
        object_significatives_bonferroni = self.selection_object_with_adjusted_error_rate(error_rate_bonferroni, df)
        significative_objects['Bonferroni'] = object_significatives_bonferroni

        object_significatives_holm = self.selection_object_with_adjusted_pvalue("Holm", df)
        significative_objects['Holm'] = object_significatives_holm

        object_significatives_sgof = self.selection_object_with_sgof("SGoF", df)
        significative_objects['SGoF'] = object_significatives_sgof

        object_significatives_benjamini_hochberg = self.selection_object_with_adjusted_pvalue("BenjaminiHochberg", df)
        significative_objects['BenjaminiHochberg'] = object_significatives_benjamini_hochberg

        return df, significative_objects

    def writing_output(self, df, significative_objects, over_or_underrepresentation, approximation_yes_or_no, yes_answers):
        df = df.sort_values(['pValueBenjaminiHochberg'])

        if approximation_yes_or_no in yes_answers:
            self.set_output_columns(1, "CountsTotal")
            df = df[self.get_output_columns()]
        else:
            df = df[self.get_output_columns()]

        if over_or_underrepresentation == 'over':
            df.to_csv(output_directory + "pValuesOf" + self.get_object_to_analyze() + "_over.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)
        elif over_or_underrepresentation == 'under':
            df.to_csv(output_directory + "pValuesOf" + self.get_object_to_analyze() + "_under.tsv", sep= "\t", float_format = '%.6f', index = True, header = True, quoting = csv.QUOTE_NONE)

        if over_or_underrepresentation == 'over':
            csvfile = open(output_directory + "significatives" + self.get_object_to_analyze() + "_over.tsv", "w", newline = "")
        elif over_or_underrepresentation == 'under':
            csvfile = open(output_directory + "significatives" + self.get_object_to_analyze() + "_under.tsv", "w", newline = "")

        writer = csv.writer(csvfile, delimiter="\t")
        writer.writerow([self.get_object_to_analyze() + 'Sidak', self.get_object_to_analyze() + 'Bonferroni', self.get_object_to_analyze() + 'Holm', self.get_object_to_analyze() + 'SGoF', self.get_object_to_analyze() + 'BenjaminiHochberg'])

        for index in range(len(significative_objects['BenjaminiHochberg'])):
            if index in range(len(significative_objects['Sidak'])):
                object_significatives_SidakValue =  significative_objects['Sidak'][index]
            else :
                object_significatives_SidakValue =  'nan'
            if index in range(len(significative_objects['Bonferroni'])):
                object_significatives_bonferroniValue =  significative_objects['Bonferroni'][index]
            else :
                object_significatives_bonferroniValue =  'nan'
            if index in range(len(significative_objects['Holm'])):
                object_significatives_holmValue =  significative_objects['Holm'][index]
            else :
                object_significatives_holmValue =  'nan'
            if index in range(len(significative_objects['SGoF'])):
                object_significatives_sgofValue =  significative_objects['SGoF'][index]
            else :
                object_significatives_sgofValue =  'nan'

            writer.writerow([object_significatives_SidakValue, object_significatives_bonferroniValue, object_significatives_holmValue, \
            object_significatives_sgofValue, significative_objects['BenjaminiHochberg'][index]])

        csvfile.close()

    def correction_bonferroni(self, df):
        pvalue_correction_bonferroni = lambda x: x * len(df.index)
        df['pValueBonferroni'] = df[self.get_statistic_method()].apply(pvalue_correction_bonferroni)

        return df

    def correction_benjamini_hochberg(self, df):
        for analyzed_object, row in df.iterrows():
            pvalue_correction_benjamini_hochberg = row[self.get_statistic_method()] * (len(df.index)/(df.index.get_loc(analyzed_object)+1))
            df.set_value(analyzed_object, 'pValueBenjaminiHochberg', pvalue_correction_benjamini_hochberg)

        return df

    def correction_holm(self, df):
        for analyzed_object, row in df.iterrows():
            pvalue_correction_holm = row[self.get_statistic_method()] * (len(df.index) - df.index.get_loc(analyzed_object))
            df.set_value(analyzed_object, 'pValueHolm', pvalue_correction_holm)

        return df

    def correction_sgof(self, df):
        F = len(df.index) * self.get_alpha()
        df = df.sort_values("pvalue_hypergeometric")
        R = 0

        R = (df['pvalue_hypergeometric'] < 0.05).sum()

        df = df.reset_index()

        object_significatives = []
        row_number = 0
        while stats.binom_test(R, len(df), p = self.get_alpha()) < self.get_alpha() and R > 0:
            df.set_value(row_number, 'pValueSGoF', 'significant')
            object_significatives.append(df.iloc[row_number][self.get_object_to_analyze()])
            R = R - 1
            row_number = row_number + 1

        df = df.set_index(self.get_object_to_analyze())

        for analyzed_object, row in df.iterrows():
            try:
                if analyzed_object not in object_significatives:
                    df.set_value(analyzed_object, 'pValueSGoF', 'nonSignificant')
            except:
                df.set_value(analyzed_object, 'pValueSGoF', 'nonSignificant')

        return df

    def error_rate_adjustement_bonferroni(self, df):
        error_rate_adjusted = self.get_alpha() / len(df.index)

        return error_rate_adjusted

    def error_rate_adjustement_sidak(self, df):
        error_rate_adjusted = (1- math.pow((1-self.get_alpha()), (1 / len(df.index))))

        return error_rate_adjusted

    def selection_object_with_adjusted_error_rate(self, error_rate, df):
        object_significatives = []

        for analyzed_object, row in df.iterrows():
            if row[self.get_statistic_method()] < error_rate :
                object_significatives.append(analyzed_object)

        return object_significatives

    def selection_object_with_adjusted_pvalue(self, method_name, df):
        object_significatives = []

        for analyzed_object, row in df.iterrows():
            if row['pValue' + method_name] < self.get_alpha() :
                object_significatives.append(analyzed_object)

        return object_significatives

    def selection_object_with_sgof(self, method_name, df):
        object_significatives = []

        for analyzed_object, row in df.iterrows():
            if row['pValue' + method_name] == 'significant' :
                object_significatives.append(analyzed_object)

        return object_significatives

    def enrichment_analysis(self):

        counts_df = pa.read_csv(temporary_directory + self.get_file_of_interest() + ".tsv", sep = "\t")
        counts_df_reference = pa.read_csv(temporary_directory + self.get_file_of_reference() + ".tsv", sep = "\t")

        counts_df = counts_df.set_index(self.get_object_to_analyze())
        counts_df_reference = counts_df_reference.set_index(self.get_object_to_analyze())

        df_joined = counts_df.join(counts_df_reference)

        yes_answers = ['yes', 'y', 'oui', 'o']
        yes_or_no = input("Is this an approximation of the reference? ")

        for analyzed_object, row in df_joined.iterrows():
            df_joined.set_value(analyzed_object, 'Percentage' + self.get_object_to_analyze() + 'InInterest', self.percentage_calculation(row['Counts'], self.get_number_of_analyzed_object_of_interest()))

        if yes_or_no in yes_answers:
            df_joined_approximation = self.counting_approximation(df_joined)
            for analyzed_object, row in df_joined_approximation.iterrows():
                df_joined_approximation.set_value(analyzed_object, 'Percentage' + self.get_object_to_analyze() + 'InReference', self.percentage_calculation(row['CountsTotal'], self.get_number_of_analyzed_object_of_reference()))

            df_joined_overrepresentation = self.hypergeometric_test_on_dataframe(df_joined_approximation, "over", 'CountsTotal')
            df_joined_overrepresentation, significative_objects = self.multiple_testing_correction(df_joined_overrepresentation)
            self.writing_output(df_joined_overrepresentation, significative_objects, "over", yes_or_no, yes_answers)

            df_joined_underrepresentation = self.hypergeometric_test_on_dataframe(df_joined_approximation, "under", 'CountsTotal')
            df_joined_underrepresentation, significative_objects = self.multiple_testing_correction(df_joined_underrepresentation)
            self.writing_output(df_joined_underrepresentation, significative_objects, "under", yes_or_no, yes_answers)

        else:
            for analyzed_object, row in df_joined.iterrows():
                df_joined.set_value(analyzed_object, 'Percentage' + self.get_object_to_analyze() + 'InReference', self.percentage_calculation(row['CountsReference'], self.get_number_of_analyzed_object_of_reference()))

            df_joined_overrepresentation = self.hypergeometric_test_on_dataframe(df_joined, "over", 'CountsReference')
            df_joined_overrepresentation, significative_objects = self.multiple_testing_correction(df_joined_overrepresentation)
            self.writing_output(df_joined_overrepresentation, significative_objects, "over", yes_or_no, yes_answers)

            df_joined_underrepresentation = self.hypergeometric_test_on_dataframe(df_joined, "under", 'CountsReference')
            df_joined_underrepresentation, significative_objects = self.multiple_testing_correction(df_joined_underrepresentation)
            self.writing_output(df_joined_underrepresentation, significative_objects, "under", yes_or_no, yes_answers)


class GOEnrichmentAnalysis(EnrichmentAnalysis):

    def __init__(self, column_name, d_go_label_to_number):
        EnrichmentAnalysis.__init__(self, column_name)
        self.output_columns.append("GOLabel")
        self.gos_labels_to_numbers = d_go_label_to_number

    def get_gos_labels_to_numbers(self):
        return self.gos_labels_to_numbers

    def tranlsation_go_number_to_go_label(self, go_numbers, d_go_label_to_number):
        go_labels = []

        for go_number in go_numbers:
            if go_number in d_go_label_to_number:
                go_labels.append(d_go_label_to_number[go_number])

        return go_labels

    def multiple_testing_correction(self, df):
        df = df.sort_values([self.get_statistic_method()])

        df = self.correction_bonferroni(df)
        df = self.correction_benjamini_hochberg(df)
        df = self.correction_holm(df)
        df = self.correction_sgof(df)

        significative_objects = {}
        translation_gos_labels_to_numbers = self.get_gos_labels_to_numbers()

        error_rate_sidak = self.error_rate_adjustement_sidak(df)
        object_significatives_Sidak = self.selection_object_with_adjusted_error_rate(error_rate_sidak, df)
        go_label_significatives_sidak = self.tranlsation_go_number_to_go_label(object_significatives_Sidak, translation_gos_labels_to_numbers)
        significative_objects['Sidak'] = go_label_significatives_sidak

        error_rate_bonferroni = self.error_rate_adjustement_bonferroni(df)
        object_significatives_bonferroni = self.selection_object_with_adjusted_error_rate(error_rate_bonferroni, df)
        go_label_significatives_bonferroni= self.tranlsation_go_number_to_go_label(object_significatives_bonferroni, translation_gos_labels_to_numbers)
        significative_objects['Bonferroni'] = go_label_significatives_bonferroni

        object_significatives_holm = self.selection_object_with_adjusted_pvalue("Holm", df)
        go_label_significatives_holm = self.tranlsation_go_number_to_go_label(object_significatives_holm, translation_gos_labels_to_numbers)
        significative_objects['Holm'] = go_label_significatives_holm

        object_significatives_sgof = self.selection_object_with_sgof("SGoF", df)
        go_label_significatives_sgof = self.tranlsation_go_number_to_go_label(object_significatives_sgof, translation_gos_labels_to_numbers)
        significative_objects['SGoF'] = go_label_significatives_sgof

        object_significatives_benjamini_hochberg = self.selection_object_with_adjusted_pvalue("BenjaminiHochberg", df)
        go_label_significatives_benjamini_hochberg = self.tranlsation_go_number_to_go_label(object_significatives_benjamini_hochberg, translation_gos_labels_to_numbers)
        significative_objects['BenjaminiHochberg'] = go_label_significatives_benjamini_hochberg

        for go, row in df.iterrows():
            if go in translation_gos_labels_to_numbers:
                df.set_value(go, 'GOLabel', translation_gos_labels_to_numbers[go])

        return df, significative_objects

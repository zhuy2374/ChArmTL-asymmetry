#
# @Author  : zhuy
# @Date    : 2026/1/10 22:05
# @Desc    : Parent-offspring similarity calculator
#
import pandas as pd
from pandas import DataFrame
from src.modules.simillarity_calculator import SimilarytiCalculator


class OffspringParentSimilarityCalculator:
    @staticmethod
    def compute_similarity(offspring_df: DataFrame, mo_df: DataFrame, fa_df: DataFrame):
        """
        Compute parent-offspring similarity.
        :param offspring_df: Offspring data.
        :param mo_df: Mother data.
        :param fa_df: Father data.
        :return:
        """

        # Build the result DataFrame.
        result_df = DataFrame(columns=['number', 'chr', 'offspring_allele_id',
                                       'offspring_TL_p75', 'offspring_reads_numer', 'offspring_tvr_len',
                                       'offspring_tvr_consensus', 'parent','parent_allele_id','parent_TL_p75',
                                       'parent_reads_numer', 'parent_tvr_len', 'parent_tvr_consensus',
                                       'similarity'
                                       ])
        # Counter for the number field.
        counter = 1

        # Result rows.
        rows = []

        # Iterate through each offspring row.
        for _, offspring in offspring_df.iterrows():
            chr = offspring['#chr']
            # Find all rows in mo_df where #chr matches chr.
            mo_rows = mo_df[mo_df['#chr'] == chr]
            # Find all rows in fa_df where #chr matches chr.
            fa_rows = fa_df[fa_df['#chr'] == chr]

            # Compute similarity with maternal telomeres.
            for _, mo in mo_rows.iterrows():
                # Compute similarity and keep 4 decimal places.
                similarity = SimilarytiCalculator.calculate_similarity(offspring['tvr_consensus'], mo['tvr_consensus'], 4)
                # Add data to rows.
                mo_data = {
                    'number': counter,
                    'chr': chr,
                    'offspring_allele_id': offspring['allele_id'],
                    'offspring_TL_p75': offspring['TL_p75'],
                    'offspring_reads_numer': offspring['reads_number'],
                    'offspring_tvr_len': offspring['tvr_len'],
                    'offspring_tvr_consensus': offspring['tvr_consensus'],
                    'parent': 'mo',
                    'parent_allele_id': mo['allele_id'],
                    'parent_TL_p75': mo['TL_p75'],
                    'parent_reads_numer': mo['reads_number'],
                    'parent_tvr_len': mo['tvr_len'],
                    'parent_tvr_consensus': mo['tvr_consensus'],
                    'similarity': similarity,
                }
                rows.append(pd.DataFrame([mo_data]))
                counter += 1

            # Compute similarity with paternal telomeres.
            for _, fa in fa_rows.iterrows():
                similarity = SimilarytiCalculator.calculate_similarity(offspring['tvr_consensus'], fa['tvr_consensus'], 4)
                fa_data = {
                    'number': counter,
                    'chr': chr,
                    'offspring_allele_id': offspring['allele_id'],
                    'offspring_TL_p75': offspring['TL_p75'],
                    'offspring_reads_numer': offspring['reads_number'],
                    'offspring_tvr_len': offspring['tvr_len'],
                    'offspring_tvr_consensus': offspring['tvr_consensus'],
                    'parent': 'fa',
                    'parent_allele_id': fa['allele_id'],
                    'parent_TL_p75': fa['TL_p75'],
                    'parent_reads_numer': fa['reads_number'],
                    'parent_tvr_len': fa['tvr_len'],
                    'parent_tvr_consensus': fa['tvr_consensus'],
                    'similarity': similarity,
                }
                rows.append(pd.DataFrame([fa_data]))
                counter += 1

        if rows:
            result_df = pd.concat(rows, ignore_index=True)
        else:
            result_df = pd.DataFrame(columns=result_df.columns)

        return result_df

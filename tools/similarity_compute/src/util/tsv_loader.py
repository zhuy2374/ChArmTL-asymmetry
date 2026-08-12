#
# @Author  : zhuy
# @Date    : 2026/1/10 17:23
# @Desc    : TSV loader utility class
#

import pandas as pd
import logging

logging.basicConfig(level=logging.INFO)


class TSVLoader:
    """
    Utility class for loading TSV files.
    """

    @staticmethod
    def load(file_path, row_filters=None, col_filters=None, logic='AND', **kwargs):
        """
        :param file_path: File path.
        :param logic: How row filters are combined.
        :param row_filters: List of row filters that return masks.
        :param col_filters: List of column operations that return processed DataFrames.
        """
        kwargs.setdefault('sep', '\t')
        df = pd.read_csv(file_path, **kwargs)

        # A. Apply row filters first.
        if row_filters:
            if logic.upper() == 'AND':
                mask = pd.Series([True] * len(df), index=df.index)
                for f in row_filters: mask &= f(df)
            else:
                mask = pd.Series([False] * len(df), index=df.index)
                for f in row_filters: mask |= f(df)
            df = df[mask]

        # B. Apply column selection or exclusion.
        if col_filters:
            for f in col_filters:
                df = f(df)

        return df


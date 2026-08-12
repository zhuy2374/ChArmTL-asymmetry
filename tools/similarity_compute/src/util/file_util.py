#
# @Author  : zhuy
# @Date    : 2026/1/10 21:55
# @Desc    : File utility class
#
import os
from pathlib import Path

import pandas as pd
from pandas import DataFrame

from src.util.tsv_loader import TSVLoader

import logging


class FileUtil:

    @staticmethod
    def load_csv(path):
        """
        Load a CSV file.
        :param path: CSV file path.
        :return: DataFrame.
        """
        # 1. Check whether the path exists.
        if not os.path.exists(path):
            raise RuntimeError(f"File not found: {path}")
        try:
            return pd.read_csv(path)
        except Exception as e:
            raise RuntimeError(f"Error loading {path}: {e}")

    @staticmethod
    def load_tsv(path):
        # 1. Check whether the path exists.
        if not os.path.exists(path):
            raise RuntimeError(f"File not found: {path}")
        try:
            return TSVLoader.load(path)
        except Exception as e:
            raise RuntimeError(f"Error loading {path}: {e}")

    @staticmethod
    def load_ped_file(path):
        """
        Load a PED-format file and validate its format.
        :param path: PED file path.
        :return: Validated DataFrame.
        """
        if not os.path.exists(path):
            raise RuntimeError(f"File not found: {path}")
        # Load the file.
        return pd.read_csv(path, sep=r'\s+', dtype=str)


    @staticmethod
    def get_real_case_path(path_str):
        """
        Return the actual filesystem path for a path string that may use incorrect casing.
        Raise FileNotFoundError if the path does not exist.
        """
        path_obj = Path(path_str)
        parts = path_obj.parts

        # Determine whether traversal starts from the root path or the current path.
        if path_obj.is_absolute():
            # On Windows, parts[0] is the drive letter, while on Linux it is '/'.
            # Windows drive letters are usually case-insensitive, so use it as the starting point.
            current_path = Path(parts[0])
            parts_to_traverse = parts[1:]
        else:
            current_path = Path(".")
            parts_to_traverse = parts

        # Traverse path components level by level.
        for part in parts_to_traverse:
            # List all entries under current_path and find one whose lowercase
            # name matches the lowercase form of part.
            found = False

            # Guard against non-directory path segments that cannot be traversed.
            if not current_path.is_dir():
                raise FileNotFoundError(f"Path segment '{current_path}' is not a directory")

            for child in current_path.iterdir():
                if child.name.lower() == part.lower():
                    current_path = child
                    found = True
                    break

            if not found:
                raise FileNotFoundError(f"Could not find '{part}' under '{current_path}'")

        return current_path


    @staticmethod
    def get_family_dfs(family_code, root_dir, row_filters, col_filters, with_5mc):
        """
        Get DataFrames for all samples in a family.
        :param family_code: Family code.
        :param root_dir: Family folder root directory.
        :param row_filters: Row filters.
        :param col_filters: Column filters.
        :param with_5mc: Whether to compute 5mc results.
        :return: s1_df, p1_df, fa_df, mo_df
        """

        # Get TSV file paths.
        if with_5mc:
            s1_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_s1/{family_code}_s1_5mc_telo_out/tlens_by_allele.tsv")
            p1_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_p1/{family_code}_p1_5mc_telo_out/tlens_by_allele.tsv")
            fa_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_fa/{family_code}_fa_5mc_telo_out/tlens_by_allele.tsv")
            mo_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_mo/{family_code}_mo_5mc_telo_out/tlens_by_allele.tsv")
        else:
            s1_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_s1/{family_code}_s1_telo_out/tlens_by_allele.tsv")
            p1_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_p1/{family_code}_p1_telo_out/tlens_by_allele.tsv")
            fa_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_fa/{family_code}_fa_telo_out/tlens_by_allele.tsv")
            mo_tsv_path = FileUtil.get_real_case_path(f"{root_dir}/{family_code}_mo/{family_code}_mo_telo_out/tlens_by_allele.tsv")

        # Read TSV files into DataFrames with configurable filters.
        s1_df = TSVLoader.load(s1_tsv_path, row_filters, col_filters)
        p1_df = TSVLoader.load(p1_tsv_path, row_filters, col_filters)
        fa_df = TSVLoader.load(fa_tsv_path, row_filters, col_filters)
        mo_df = TSVLoader.load(mo_tsv_path, row_filters, col_filters)

        return s1_df, p1_df, fa_df, mo_df

    @staticmethod
    def write_dataframe_to_csv(df: DataFrame, path, sort: bool = True):
        if sort:
            # Sort by similarity in descending order.
            df = FileUtil.sort_dataframe(df, 'offspring_allele_id', 'similarity', index="number")
        df.to_csv(path, index=False, encoding="utf-8")

    @staticmethod
    def write_dataframe_to_tsv(df: DataFrame, path, sort: bool = True):
        if sort:
            # Sort by similarity in descending order.
            df = FileUtil.sort_dataframe(df, 'offspring_allele_id', 'similarity', index="number")
        df.to_csv(path, index=False, encoding="utf-8", sep="\t")


    @staticmethod
    def sort_dataframe(df: DataFrame, group_columns: str, sort_columns:  str, ascending: bool = False, index: str = None) -> DataFrame:
        """
        Sort a DataFrame.
        :param df: Data to sort.
        :param group_columns: Grouping column.
        :param sort_columns: Sort column.
        :param ascending: Whether to sort in ascending order.
        :param index: Sequence number column.
        :return: Sorted DataFrame.
        """
        df[group_columns] = pd.Categorical(df[group_columns], categories=pd.unique(df[group_columns]), ordered=True)
        df = df.sort_values([group_columns, sort_columns], ascending=[True, ascending])

        # Reset sequence numbers if requested.
        if index is not None:
            df[index] = range(1, len(df) + 1)

        # Reset the index.
        df = df.reset_index(drop=True)

        return df

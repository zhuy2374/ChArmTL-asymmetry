#
# @Author  : zhuy
# @Date    : 2026/1/10 20:07
# @Desc    : Similarity calculation script
#
import argparse
from pathlib import Path
import logging
import sys
import os

from tqdm import tqdm
# Get the project root path and insert it at the start of sys.path so project modules can be imported.
project_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, project_root)

from src.modules.family_mapper import FamilyMapper
from src.modules.offspring_parent_similarity_calculator import OffspringParentSimilarityCalculator
from src.util.file_util import FileUtil
from src.util.tsv_filters import column_not_contains, drop_columns, fast_assign, column_not_nan

logging.basicConfig(level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description="Collect family directories named <familyID>_<individual>.")
    parser.add_argument("root_dir", help="root directory containing Telogator2 results for samples", type=str)
    parser.add_argument("--with_5mc", action="store_true", help="specify if the Telogator2 results were generated from 5mc BAM files")
    parser.add_argument("--save_csv", action="store_true", help="save the output to a CSV file",default= True)
    parser.add_argument("--save_tsv", action="store_true", help="save the output to a TSV file")
    parser.add_argument("--sort", action="store_true", help="enable sorting of the results")
    args = parser.parse_args()

    try:

        base_dir = args.root_dir

        logging.info(f"Starting to process root directory: {base_dir}")
        # Count families under the root directory.
        family_dict = FamilyMapper.collect_family_dirs(root=Path(base_dir), return_absolute=True)
        logging.info(f"Directory scan complete. Found {len(family_dict)} families.")

        # Define filters.
        row_filters = [
            # Filter rows.
            column_not_contains('#chr', ','),
            # Filter out rows where tvr_consensus is NaN.
            column_not_nan('tvr_consensus'),
        ]
        col_filters = [
            # Drop unused columns.
            drop_columns(["position", "ref_samp", "supporting_reads"]),
            # Add dynamic columns.
            fast_assign(
                reads_number=lambda x: x['read_TLs'].str.split(",").str.len()
            )
        ]
        logging.info("Starting to process family data...")

        # Iterate through each family.
        for key, value in tqdm(family_dict.items(), total=len(family_dict), desc="Processing family data"):
            # Get DataFrames.
            s1_df, p1_df, fa_df, mo_df = FileUtil.get_family_dfs(key, base_dir, row_filters, col_filters, with_5mc=args.with_5mc)

            # Compute similarity between s1 and parents.
            s1_similarity_result_df = OffspringParentSimilarityCalculator.compute_similarity(s1_df, mo_df, fa_df)

            # Compute similarity between p1 and parents.
            p1_similarity_result_df = OffspringParentSimilarityCalculator.compute_similarity(p1_df, mo_df, fa_df)

            if args.save_csv:
                if args.with_5mc:
                    FileUtil.write_dataframe_to_csv(s1_similarity_result_df, Path(f"{base_dir}/{key}_s1/{key}_s1_5mc_similarity.csv"), sort=args.sort)
                    FileUtil.write_dataframe_to_csv(p1_similarity_result_df, Path(f"{base_dir}/{key}_p1/{key}_p1_5mc_similarity.csv"), sort=args.sort)
                else:
                    FileUtil.write_dataframe_to_csv(s1_similarity_result_df, Path(f"{base_dir}/{key}_s1/{key}_s1_similarity.csv"), sort=args.sort)
                    FileUtil.write_dataframe_to_csv(p1_similarity_result_df, Path(f"{base_dir}/{key}_p1/{key}_p1_similarity.csv"), sort=args.sort)

            if args.save_tsv:
                if args.with_5mc:
                    FileUtil.write_dataframe_to_tsv(s1_similarity_result_df, Path(f"{base_dir}/{key}_s1/{key}_s1_5mc_similarity.tsv"), sort=args.sort)
                    FileUtil.write_dataframe_to_tsv(p1_similarity_result_df, Path(f"{base_dir}/{key}_p1/{key}_p1_5mc_similarity.tsv"), sort=args.sort)
                else:
                    FileUtil.write_dataframe_to_tsv(s1_similarity_result_df, Path(f"{base_dir}/{key}_s1/{key}_s1_similarity.tsv"), sort=args.sort)
                    FileUtil.write_dataframe_to_tsv(p1_similarity_result_df, Path(f"{base_dir}/{key}_p1/{key}_p1_similarity.tsv"), sort=args.sort)

        logging.info("Processing complete.")
        if args.save_csv:
            logging.info("All results have been saved to CSV.")

        if args.save_tsv:
            logging.info("All results have been saved to TSV.")

    except Exception as e:
        logging.error(f"Error while processing data: {e}")
        return


if __name__ == "__main__":
    main()

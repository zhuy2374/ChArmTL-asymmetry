#
# @Author  : zhuy
# @Date    : 2026/1/13 17:44
# @Desc    : PED-based family similarity calculation script
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

from src.util.tsv_loader import TSVLoader
from src.modules.offspring_parent_similarity_calculator import OffspringParentSimilarityCalculator
from src.util.file_util import FileUtil
from src.util.tsv_filters import column_not_contains, drop_columns, fast_assign, column_not_nan

logging.basicConfig(level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description="Collect family directories named <familyID>_<individual>.")
    parser.add_argument("root_dir", help="root directory containing Telogator2 results for samples", type=str)
    parser.add_argument("--family_list", help="file containing family codes", type=str)
    parser.add_argument("--with_5mc", action="store_true", help="specify if the Telogator2 results were generated from 5mc BAM files")
    parser.add_argument("--save_csv", action="store_true", help="save the output to a CSV file",default= True)
    parser.add_argument("--save_tsv", action="store_true", help="save the output to a TSV file")
    parser.add_argument("--sort", action="store_true", help="enable sorting of the results")
    args = parser.parse_args()

    try:
        base_dir = args.root_dir

        logging.info(f"Starting to process root directory: {base_dir}")

        # Parse the PED-format file.
        ped_df = FileUtil.load_ped_file(args.family_list)

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
        # Iterate through each individual.
        for _, row in tqdm(ped_df.iterrows(), total=len(ped_df), desc="Processing family data"):
            offspring = row['IID']
            mo = row['MID']
            fa = row['PID']
            # Skip parents.
            if str(mo) == "0" and str(fa) == "0":
                continue

            offspring_path = f"{base_dir}/{offspring}/{offspring}{'_5mc_telo_out' if args.with_5mc else '_telo_out'}/tlens_by_allele.tsv"
            mo_path = f"{base_dir}/{mo}/{mo}{'_5mc_telo_out' if args.with_5mc else '_telo_out'}/tlens_by_allele.tsv"
            fa_path = f"{base_dir}/{fa}/{fa}{'_5mc_telo_out' if args.with_5mc else '_telo_out'}/tlens_by_allele.tsv"

            offspring_df = TSVLoader.load(Path(offspring_path), row_filters, col_filters)
            mo_df = TSVLoader.load(Path(mo_path), row_filters, col_filters)
            fa_df = TSVLoader.load(Path(fa_path), row_filters, col_filters)

            result_df = OffspringParentSimilarityCalculator.compute_similarity(offspring_df, mo_df, fa_df)

            if args.save_csv:
                if args.with_5mc:
                    FileUtil.write_dataframe_to_csv(result_df, Path(f"{base_dir}/{offspring}/{offspring}_5mc_similarity.csv"), sort=args.sort)
                else:
                    FileUtil.write_dataframe_to_csv(result_df, Path(f"{base_dir}/{offspring}/{offspring}_similarity.csv"), sort=args.sort)

            if args.save_tsv:
                if args.with_5mc:
                    FileUtil.write_dataframe_to_tsv(result_df, Path(f"{base_dir}/{offspring}/{offspring}_5mc_similarity.tsv"), sort=args.sort)
                else:
                    FileUtil.write_dataframe_to_tsv(result_df, Path(f"{base_dir}/{offspring}/{offspring}_similarity.tsv"), sort=args.sort)

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

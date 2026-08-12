### Pedigree-Based Computational Tool for Telomere Variant Repeat (TVR) Sequence Similarity

#### Getting Started

1. Install **uv** and **python 3.12**:

   ```shell
   # Install the uv package manager
   curl -LsSf https://astral.sh/uv/install.sh | sh
   
   # Install python 3.12
   uv python install 3.12
   ```

2. Sync dependencies:

   ```shell
   uv sync
   ```

3. Run the script using **uv**:

   ```shell
   # If individual Telogator2 result files are named as "<family_id>_<role>" (e.g., "FJ0018_fa"), run:
   uv run ./scripts/run_compute_similarity.py "<root_dir>"
   
   # If individual Telogator2 result files do not follow the family naming convention, run:
   uv run ./scripts/run_compute_similarity_use_ped.py "<root_dir>"
   ```

   #### Parameter Description：

   **run_compute_similarity.py**：
   
   |  Parameter     |                             Description                             |
   |:--------------:|:-------------------------------------------------------------------:|
   |    root_dir    |      Root directory containing Telogator2 results for samples       |
   |   --with_5mc   | Specify if the Telogator2 results were generated from 5mc BAM files |
   |   --save_csv   |                    Save the output to a CSV file                    |
   |   --save_tsv   |                    Save the output to a TSV file                    |
   |     --sort     |                    Enable sorting of the results                    |
   
   **run_compute_similarity_use_ped.py**：
   
   |   Parameter   |                             Description                             |
   | :-----------: |:-------------------------------------------------------------------:|
   |   root_dir    |     Root directory containing Telogator2 results for samples        |
   |  --with_5mc   | Specify if the Telogator2 results were generated from 5mc BAM files |
   |  --save_csv   |                    Save the output to a CSV file                    |
   |  --save_tsv   |                    Save the output to a TSV file                    |
   | --family_list |               Path to the simple PED (pedigree) file                |
   |    --sort     |                    Enable sorting of the results                    |
   
   #### Simple PED File Format
   
   If you are using the `run_compute_similarity_use_ped.py` script, you must provide a tab-separated text file via the `--family_list` parameter.
   
   **Column Definitions:**
   * `FID`: Family ID
   * `IID`: Individual ID
   * `PID`: Paternal ID (`0` if unknown)
   * `MID`: Maternal ID (`0` if unknown)

   **An example (`family_list.ped`):**
   
   ```text
   FID	IID	PID	MID
   FJ0018	FJ0018_fa	0	0
   FJ0018	FJ0018_mo	0	0
   FJ0018	FJ0018_p1	FJ0018_fa	FJ0018_mo
   FJ0018	FJ0018_s1	FJ0018_fa	FJ0018_mo
   ```

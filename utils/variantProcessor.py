"""
VariantProcessor - A utility class for processing variant presence/absence data

This module filters variants based on their frequency across samples.
"""

import os
import logging
import pandas as pd
from typing import Optional, Dict

# Set up logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class VariantProcessor:
    """
    A class to filter variants based on frequency threshold.
    """

    def __init__(self, input_dir: str = '00_input_files'):
        self.input_dir = input_dir
        self.raw_df = None
        self.filtered_df = None

    def load_data(self, filename: str = 'Effector_variants_PAV_output.txt') -> pd.DataFrame:
        """
        Load the variant presence/absence matrix from file.

        Args:
            filename: File name to load (relative to input_dir or full path)

        Returns:
            Loaded DataFrame
        """
        file_path = filename if os.path.isabs(filename) else os.path.join(self.input_dir, filename)
        try:
            self.raw_df = pd.read_csv(file_path, sep='\t')
            logger.info(f"Loaded variant data from {file_path} with shape {self.raw_df.shape}")
            return self.raw_df
        except Exception as e:
            logger.error(f"Failed to load variant data from {file_path}: {e}")
            raise

    def filter_by_variant_frequency(self, min_var: int = 5) -> pd.DataFrame:
        """
        Filter variants that occur in fewer samples than the minimum frequency.

        Args:
            min_var: Minimum frequency threshold to retain a variant

        Returns:
            Filtered DataFrame
        """
        if self.raw_df is None:
            logger.error("Raw data not loaded. Call load_data() before filtering.")
            raise ValueError("Data must be loaded before filtering.")

        try:
            variant_sums = self.raw_df.iloc[:, 1:].sum()
            filtered_variants = variant_sums[variant_sums > min_var]
            self.filtered_df = self.raw_df[['ID'] + list(filtered_variants.index)].copy()

            logger.info(f"{len(filtered_variants)} variants retained (min_var > {min_var})")
            return self.filtered_df
        except Exception as e:
            logger.error(f"Failed to filter variants by frequency: {e}")
            raise

    def save_processed_data(self, output_file: str = '01_intermediate_files/0_filtered_combined_variants.txt') -> None:
        """
        Save the processed DataFrame to a file.

        Args:
            output_file: Path to output file
        """
        if self.filtered_df is None:
            logger.error("No processed data available to save.")
            raise ValueError("No processed data to save.")

        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        self.filtered_df.to_csv(output_file, sep='\t', index=False)
        logger.info(f"Processed variant data saved to: {output_file}")

    def process_data(self, var_data: pd.DataFrame) -> Dict[str, pd.DataFrame]:
        """
        Placeholder for advanced processing of variant data.

        Args:
            var_data: A DataFrame of binary presence/absence calls for each variant

        Returns:
            Dictionary mapping output names to processed DataFrames
        """
        try:
            # Future extension: e.g., clustering, dimensionality reduction, etc.
            result = {"processed_variants": var_data.copy()}
            logger.info(f"(placeholder): processed {len(var_data.columns) - 1} variants for {len(var_data)} samples.")
            return result
        except Exception as e:
            logger.error(f"Data processing failed: {e}")
            raise

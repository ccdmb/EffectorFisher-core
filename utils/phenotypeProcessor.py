"""
PhenotypeProcessor - A utility class for processing phenotype data

This module provides functionality for processing phenotype data,
either quantitative or qualitative, and transforming it as needed.
"""

import os
import logging
import pandas as pd
from typing import Dict, Optional

# Set up logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)  # Change to DEBUG for more verbosity
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class PhenotypeProcessor:
    """
    A class to process phenotype data from various formats.
    """

    def __init__(self, input_dir: str = '00_input_files'):
        self.input_dir = input_dir
        self.quantitative_data = None
        self.qualitative_data = None
        self.processed_traits = {}

    def load_quantitative_data(self, filename: str = '0_phenotype_data_quantitative.txt') -> pd.DataFrame:
        file_path = os.path.join(self.input_dir, filename)
        try:
            df = pd.read_csv(file_path, sep='\t', index_col=0)
            df.reset_index(inplace=True)
            df.rename(columns={df.columns[0]: 'ID'}, inplace=True)
            self.quantitative_data = df
            return df
        except Exception as e:
            logger.error(f"Failed to load quantitative data from {file_path}: {e}")
            raise

    def load_qualitative_data(self, filename: str = '0_phenotype_data_qualitative.txt') -> pd.DataFrame:
        file_path = os.path.join(self.input_dir, filename)
        try:
            df = pd.read_csv(file_path, sep='\t')
            if df.columns[0] != 'ID':
                df.columns = ['ID'] + df.columns[1:].tolist()
            self.qualitative_data = df
            return df
        except Exception as e:
            logger.error(f"Failed to load qualitative data from {file_path}: {e}")
            raise

    def convert_quantitative_to_qualitative(self) -> pd.DataFrame:
        if self.quantitative_data is None:
            logger.error("Quantitative data not loaded. Call load_quantitative_data first.")
            raise ValueError("Quantitative data not loaded. Call load_quantitative_data first.")

        try:
            converted_df = self.quantitative_data.copy()
            for column in converted_df.columns[1:]:
                median_value = converted_df[column].median()
                converted_df[column] = converted_df[column].apply(
                    lambda x: 'low' if x < median_value else 'high'
                )
            self.qualitative_data = converted_df
            return converted_df
        except Exception as e:
            logger.error(f"Failed to convert quantitative data to qualitative: {e}")
            raise

    def split_traits_to_separate_dataframes(self) -> Dict[str, pd.DataFrame]:
        if self.qualitative_data is None:
            logger.error("Qualitative data not loaded. Call load_qualitative_data or convert_quantitative_to_qualitative first.")
            raise ValueError("Qualitative data not loaded. Call load_qualitative_data or convert_quantitative_to_qualitative first.")

        try:
            result = {}
            for i, col_name in enumerate(self.qualitative_data.columns[1:], start=1):
                subset_df = self.qualitative_data[['ID', col_name]].copy()
                subset_df.columns = ['ID', 'disease']
                cleaned_df = self._clean_missing_disease_data(subset_df)
                result[f'trait_{i}'] = cleaned_df

            self.processed_traits = result
            return result
        except Exception as e:
            logger.error(f"Failed to split traits into separate DataFrames: {e}")
            raise

    def _clean_missing_disease_data(self, df: pd.DataFrame) -> pd.DataFrame:
        try:
            return df[df.iloc[:, 1].notna()].copy()
        except Exception as e:
            logger.error(f"Failed to clean missing disease data: {e}")
            raise

    def save_processed_data(self, output_dir: str = 'output') -> Dict[str, pd.DataFrame]:
        """
        Save qualitative and per-trait dataframes to disk, and return them.

        Args:
            output_dir: Output directory to save the files

        Returns:
            Dictionary mapping file paths to their corresponding DataFrames
        """
        output_data = {}
        os.makedirs(output_dir, exist_ok=True)

        try:
            # Save full qualitative dataset
            if self.qualitative_data is not None:
                qual_path = os.path.join(output_dir, '0_phenotype_data_qualitative.txt')
                self.qualitative_data.to_csv(qual_path, sep='\t', index=False)
                output_data[qual_path] = self.qualitative_data.copy()

            # Save cleaned individual trait dataframes
            for i, (trait_name, trait_df) in enumerate(self.processed_traits.items(), start=1):
                cleaned_df = self._clean_missing_disease_data(trait_df)
                output_path = os.path.join(output_dir, f'1_data{i}.txt')
                cleaned_df.to_csv(output_path, sep='\t', index=False)
                output_data[output_path] = cleaned_df

            return output_data
        except Exception as e:
            logger.error(f"Failed to save processed data: {e}")
            raise

    def process_data(self, data_type: str = 'quantitative') -> Dict[str, pd.DataFrame]:
        """
        Run the full processing pipeline: load data, convert (if quantitative), and split.

        Args:
            data_type: 'quantitative' or 'qualitative'

        Returns:
            Dictionary of processed trait DataFrames
        """
        try:
            if data_type == 'quantitative':
                self.load_quantitative_data()
                self.convert_quantitative_to_qualitative()
            else:
                self.load_qualitative_data()

            return self.split_traits_to_separate_dataframes()
        except Exception as e:
            logger.error(f"Data processing failed: {e}")
            raise

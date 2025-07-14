import os
import pandas as pd
import logging

# Logger setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class ProcessPredector:
    """
    ProcessPredector - A utility class for processing Predector input

    This module merges in-memory Predector results with output from EffectorFisher.
    """

    def __init__(self, input_dir: str = '00_input_files'):
        self.input_dir = input_dir
        self.predector_df = pd.DataFrame()
        self.fisher_df = pd.DataFrame()
        self.merged_df = pd.DataFrame()

    def load_data_predector(self):
        """Load Predector results from a default file in the input directory."""
        predector_file = os.path.join(self.input_dir, "0_predector_results.txt")
        try:
            self.predector_df = pd.read_csv(predector_file, sep="\t")
            logger.info(f"Loaded Predector file: {predector_file}")
        except Exception as e:
            logger.error(f"Failed to load Predector file: {e}")
            raise


    def merge_with_fisher(self, fisher_object):
        """
        Merge Fisher result DataFrame (already processed) with Predector results.
        """
        if not hasattr(fisher_object, "merged_with_locus_df"):
            raise ValueError("Fisher object does not have 'merged_with_locus_df'.")

        self.fisher_df = fisher_object.merged_with_locus_df.copy()
        self.fisher_df = self._move_locus_id_first(self.fisher_df)
        self.predector_df = self._move_locus_id_first(self.predector_df)

        try:
            self.merged_df = self.fisher_df.merge(self.predector_df, on="locus_id", how="left")
            self.merged_df.fillna('NA', inplace=True)
            logger.info("Step 10 completed: Predector results merged with Fisher output.")
        except Exception as e:
            logger.error(f"Error during merge: {e}")
            raise


    def _move_locus_id_first(self, df: pd.DataFrame) -> pd.DataFrame:
        """Ensure 'locus_id' column is the first column."""
        if 'locus_id' in df.columns:
            cols = ['locus_id'] + [col for col in df.columns if col != 'locus_id']
            return df[cols]
        return df

    def save_processed_data(pred, fisher, output_path: str):
        """
        Save merged data after ensuring the column order matches:
        - All fisher columns (from fisher.merged_with_locus_df)
        - Followed by remaining predector columns
        """
        if pred.merged_df.empty:
            logger.warning("No merged data to save. Please run merge_with_fisher first.")
            return

        try:
            # Reorder columns
            cols_fisher = list(fisher.merged_with_locus_df.columns)
            cols_predector = [col for col in pred.predector_df.columns if col != 'locus_id']
            pred.merged_df = pred.merged_df[cols_fisher + cols_predector]

            # Save to file
            pred.merged_df.to_csv(output_path, sep="\t", index=False)
            logger.info(f"Merged dataset with reordered columns saved to: {output_path}")
        except Exception as e:
            logger.error(f"Failed to save reordered merged data: {e}")
            raise

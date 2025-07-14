import os
import logging
import pandas as pd
from typing import Dict

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class EffectorAnnotator:
    """
    EffectorAnnotator - Annotates merged Fisher-Predector data with known effectors
    and maps cultivar names to p-value columns.
    """

    def __init__(self, fisher_predector_df: pd.DataFrame, phenotype_df: pd.DataFrame, input_dir: str = "00_input_files"):
        self.fisher_predector_df = fisher_predector_df.copy()
        self.phenotype_df = phenotype_df.copy()
        self.input_dir = input_dir
        self.annotated_df = pd.DataFrame()

    def add_known_effectors(self, known_effector_file: str = "known_effector.txt") -> pd.DataFrame:
        """
        Annotates the DataFrame with known effectors based on locus_id,
        rearranges key columns, and renames p-value columns using cultivar names.
        """
        known_effectors_dict = {}

        known_path = os.path.join(self.input_dir, known_effector_file)
        if os.path.exists(known_path):
            with open(known_path, "r") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        known_effectors_dict[parts[0]] = parts[1]

        self.fisher_predector_df['known_effector'] = self.fisher_predector_df['locus_id'].map(
            known_effectors_dict).fillna("")

        logger.info("Known effectors annotated.")

        # Run internal transformations
        self._rearrange_columns()
        self._rename_pvalue_columns_with_cultivars()

        return self.annotated_df

    def _rearrange_columns(self) -> None:
        """
        Rearranges columns to bring key values near 'effector_score'.
        """
        important_cols = ["residue_number", "aa_c_number", "effector_matches", "known_effector"]

        if "effector_score" in self.fisher_predector_df.columns:
            score_idx = self.fisher_predector_df.columns.get_loc("effector_score")
            reordered_cols = list(self.fisher_predector_df.columns)

            for col in reversed(important_cols):
                if col in reordered_cols:
                    reordered_cols.remove(col)
                    reordered_cols.insert(score_idx + 1, col)

            self.annotated_df = self.fisher_predector_df[reordered_cols].copy()
            logger.info("Columns rearranged.")
        else:
            self.annotated_df = self.fisher_predector_df.copy()
            logger.warning("effector_score column not found; skipping rearrangement.")

    def _rename_pvalue_columns_with_cultivars(self) -> None:
        """
        Replaces p-value column names with cultivar names from phenotype data.
        """
        pval_cols = [col for col in self.annotated_df.columns if col.startswith("p-value-")]
        cultivar_names = list(self.phenotype_df.columns[1:])  # skip ID

        if len(pval_cols) != len(cultivar_names):
            logger.warning("Mismatch between number of p-value columns and cultivar names. Skipping renaming.")
            return

        column_mapping = dict(zip(pval_cols, cultivar_names))
        self.annotated_df.rename(columns=column_mapping, inplace=True)
        logger.info("Cultivar names used for p-value columns.")

    def save(self, output_path: str = "02_results/complete_isoform_list.txt") -> None:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        self.annotated_df.to_csv(output_path, sep="\t", index=False)
        logger.info(f"Annotated dataset saved to: {output_path}")

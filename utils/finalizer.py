import os
import logging
import pandas as pd
from typing import Optional

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class Finalizer:
    """
    Finalizer - Performs final ranking, filtering, and known effector extraction.

    Steps:
    - Deduplicate by locus_id and rank by effector_score.
    - Apply filters on cysteine count, residue length, p-value, and score.
    - Extract and rank known effectors after filtering.
    """

    def __init__(self, annotated_df: pd.DataFrame):
        self.annotated_df = annotated_df.copy()
        self.loci_ranked_df = pd.DataFrame()
        self.filtered_df = pd.DataFrame()
        self.known_ranked_df = pd.DataFrame()

    def _finalize_loci(self) -> None:
        sorted_df = self.annotated_df.sort_values(by="p-value_lowest", ascending=True)
        dedup_df = sorted_df.drop_duplicates(subset='locus_id', keep='first')
        self.loci_ranked_df = dedup_df.sort_values(by="effector_score", ascending=False).copy()
        logger.info("Locus deduplication and ranking done.")

    def _apply_filters(self, cyst_threshold: float, max_residue_length: int,
                       min_effector_score: float, max_p_value: float) -> None:
        self.filtered_df = self.loci_ranked_df[
            (self.loci_ranked_df['aa_c_number'] >= cyst_threshold) &
            (self.loci_ranked_df['residue_number'] <= max_residue_length) &
            (self.loci_ranked_df['effector_score'] >= min_effector_score) &
            (self.loci_ranked_df['p-value_lowest'] <= max_p_value)
        ].copy()
        logger.info("Dataset filtered based on thresholds.")

    def rank_known_effectors(self, cyst_threshold: float, max_residue_length: int,
                             min_effector_score: float, max_p_value: float) -> Optional[pd.DataFrame]:
        try:
            self._finalize_loci()
            self._apply_filters(
                cyst_threshold=cyst_threshold,
                max_residue_length=max_residue_length,
                min_effector_score=min_effector_score,
                max_p_value=max_p_value
            )

            df = self.filtered_df.copy()
            df['known_effector'] = df['known_effector'].astype(str)
            df = df.sort_values(by="effector_score", ascending=False)
            df['Ranking'] = range(1, len(df) + 1)

            known = df[
                (df['known_effector'].str.strip() != "") &
                (df['known_effector'].str.lower() != 'nan')
            ].copy()

            if known.empty:
                logger.info("No known effectors found after filtering.")
                return None

            cols = ['Ranking', 'known_effector', 'locus_id', 'effector_score', 'p-value_lowest']
            self.known_ranked_df = known[cols]
            logger.info("Known effectors ranked.")
            return self.known_ranked_df
        except Exception as e:
            logger.error(f"Finalization failed: {e}")
            raise

    def save(self, output_dir: str = "02_results") -> None:
        os.makedirs(output_dir, exist_ok=True)
        self.loci_ranked_df.to_csv(os.path.join(output_dir, "complete_loci_list.txt"), sep="\t", index=False)
        self.filtered_df.to_csv(os.path.join(output_dir, "filtered_loci_list.txt"), sep="\t", index=False)

        if not self.known_ranked_df.empty:
            self.known_ranked_df.to_csv(os.path.join(output_dir, "known_effectors_ranking.txt"), sep="\t", index=False)

    def save_filtered_loci_only(self, output_dir: str = "02_results") -> None:
        """
        Save only the filtered loci list to 'filtered_loci_list.txt'
        (Used when --save is not passed but minimal output is needed)
        """
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, "filtered_loci_list.txt")

        if not self.filtered_df.empty:
            self.filtered_df.to_csv(output_path, sep="\t", index=False)
            logger.info(f"Filtered loci list saved to {output_path}")
        else:
            logger.warning("Filtered dataframe is empty. Nothing to save.")

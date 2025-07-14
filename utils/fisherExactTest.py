"""
FisherExactTest - A utility class for generating contingency and hypergeometric tables

This module merges in-memory phenotype trait data with filtered variant data,
computes contingency tables, calculates Fisher's exact test p-values,
and generates merged output for summary analysis.
"""

import os
import math
import logging
import pandas as pd
from typing import Dict

# Logger setup
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)


class FisherExactTest:

    def __init__(self, trait_data: Dict[str, pd.DataFrame], variant_df: pd.DataFrame, output_dir: str = '01_intermediate_files'):
        self.trait_data = trait_data
        self.variant_df = variant_df
        self.output_dir = output_dir
        self.hypergeo_tables: Dict[str, pd.DataFrame] = {}
        self.p_value_results: Dict[str, pd.DataFrame] = {}
        self.merged_summary_df: pd.DataFrame = pd.DataFrame()
        self.merged_with_locus_df: pd.DataFrame = pd.DataFrame()

    def generate(self) -> Dict[str, pd.DataFrame]:

        for i, (trait_name, trait_df) in enumerate(self.trait_data.items(), start=1):
            merged_df = pd.merge(trait_df, self.variant_df, on='ID', how='left')
            file_number = str(i)
            logger.info(f"Merged trait '{trait_name}' with variant data")
            hyper_df = self._create_contingency_and_hypergeo_tables(merged_df)
            self.hypergeo_tables[file_number] = hyper_df

        logger.info("Hypergeometric table generation complete.")
        return self.hypergeo_tables

    def _create_contingency_and_hypergeo_tables(self, data: pd.DataFrame) -> pd.DataFrame:

        variant_columns = data.columns[2:]
        counts = {(disease, presence): {var: 0 for var in variant_columns}
                  for disease in ['high', 'low'] for presence in [0, 1]}

        for _, row in data.iterrows():
            disease_status = row['disease'].strip().split()[0].lower() # get the type of disease, high or low

            for var in variant_columns:
                # loop over the genes
                value = row[var] # gets the presence absence, 1 or 0

                if pd.notna(value):
                    try:
                        val_int = int(value)
                        if val_int in [0, 1]:
                            counts[(disease_status, val_int)][var] += 1
                    except ValueError:
                        continue

        contingency_df = pd.DataFrame.from_dict(counts, orient='index').transpose()
        contingency_df.columns = [f"{disease}-{presence}" for (disease, presence) in contingency_df.columns]
        contingency_df = contingency_df[sorted(contingency_df.columns, key=lambda x: (x.split('-')[1], x.split('-')[0]))]

        rename_map = {'high-0': 'c', 'high-1': 'a', 'low-0': 'd', 'low-1': 'b'}
        hypergeo_df = contingency_df.rename(columns=rename_map)
        return hypergeo_df[['c', 'd', 'a', 'b']]

    def calculate_factorials(self, n: int):
        '''
        Modified factorial calculation
        '''
        fact = [0.0] * (n + 1)
        fact[0] = 1

        for i in range(2, n + 1):
            fact[i] = fact[i - 1] + math.log(i)
        return fact

    def compute_p_values(self, max_factorial: int = 100000) -> Dict[str, pd.DataFrame]:
        fact_table = self.calculate_factorials(max_factorial)

        for file_number, df in self.hypergeo_tables.items():
            records = []
            for variant, row in df.iterrows():
                try:
                    a = round(float(row['a']))
                    c = round(float(row['c']))
                    d = round(float(row['d']))
                    b = round(float(row['b']))
                    n = a + b + c + d

                    hyp1 = (
                        fact_table[a + b] +
                        fact_table[c + d] +
                        fact_table[a + c] +
                        fact_table[b + d] -
                        (fact_table[n] + fact_table[a] + fact_table[b] + fact_table[c] + fact_table[d])
                    )
                    p_value = math.exp(hyp1)

                    records.append([variant, c, d, a, b, p_value])
                except Exception as e:
                    logger.warning(f"Skipping variant '{variant}': {e}")
                    records.append([variant, 'ERROR', str(e)])

            result_df = pd.DataFrame(records, columns=['variant', 'c', 'd', 'a', 'b', 'p-value'])
            self.p_value_results[file_number] = result_df

        logger.info("Fisher exact p-values computed.")
        return self.p_value_results

    def merge_and_compute_lowest_p_value(self) -> pd.DataFrame:
        
        merged_data = pd.DataFrame()
        for file_number in sorted(self.p_value_results.keys(), key=int):
            df = self.p_value_results[file_number][['variant', 'p-value']].copy()
            df.rename(columns={'p-value': f'p-value-{file_number}'}, inplace=True)
            if merged_data.empty:
                merged_data = df
            else:
                merged_data = pd.merge(merged_data, df, on='variant', how='outer')

        p_value_cols = [col for col in merged_data.columns if col.startswith('p-value')]
        merged_data["p-value_lowest"] = merged_data[p_value_cols].min(axis=1)
        self.merged_summary_df = merged_data

        logger.info("Merged lowest p-values.")
        return merged_data

    def add_locus_id_column(self) -> pd.DataFrame:
        try:
            df = self.merged_summary_df.copy()
            df['locus_id'] = df['variant']
            df['locus_id'] = df['locus_id'].str.replace(r'(_\d+)$', '', regex=True)
            df['locus_id'] = df['locus_id'].str.replace(r'(?i)([A-D]+)$', '', regex=True)
            cols = ['locus_id'] + [col for col in df if col != 'locus_id']
            df = df[cols]
            self.merged_with_locus_df = df

            logger.info("Locus IDs extracted.")
            return df
        except Exception as e:
            logger.error(f"Failed to process locus IDs: {e}")
            raise

    def save_processed_data(self):
        os.makedirs(self.output_dir, exist_ok=True)

        for file_number, df in self.p_value_results.items():
            output_file = os.path.join(self.output_dir, f'5_output_hypergeo_data{file_number}.txt')
            df.to_csv(output_file, sep='\t', index=False)
            logger.info(f"Saved: {output_file}")

        summary_file = os.path.join(self.output_dir, '6_merged_lowest_p-value.txt')
        self.merged_summary_df.to_csv(summary_file, sep='\t', index=False)
        logger.info(f"Saved merged summary: {summary_file}")

        locus_file = os.path.join(self.output_dir, '7_merged_lowest_p-value_with_locus_id.txt')
        self.merged_with_locus_df.to_csv(locus_file, sep='\t', index=False)
        logger.info(f"Saved merged locus summary: {locus_file}")

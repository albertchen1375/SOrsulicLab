import pandas as pd
from scipy import stats
import numpy as np
from statsmodels.stats.multitest import multipletests


def perform_correlations(features, dataset):
    feature_cell_lines = features.iloc[:, 0]
    dataset_cell_lines = dataset.iloc[:, 0]
    
    common_cell_lines = set(feature_cell_lines).intersection(set(dataset_cell_lines))
    common_cell_lines = list(common_cell_lines)
    
    features = features[features.iloc[:, 0].isin(common_cell_lines)].reset_index(drop=True)
    dataset = dataset[dataset.iloc[:, 0].isin(common_cell_lines)].reset_index(drop=True)
    
    features = features.sort_values(by=features.columns[0]).reset_index(drop=True)
    dataset = dataset.sort_values(by=dataset.columns[0]).reset_index(drop=True)
    
    print(features)
    print(dataset)
    
    # Initialize a list to hold results
    results = []

    # Iterate over each column in the features DataFrame
    for col_name in features.columns:
        # Skip the column if it contains "ID" in the name
        if "depmap_id" in col_name or "nuclei_count" in col_name or "cell_line_display_name" in col_name:
            continue
        
        col_data = features[col_name].values

        # Iterate over each gene in the gene expression DataFrame
        for name in dataset.columns:
            # Skip the first column if it's an index or identifier
            if any(keyword in name for keyword in ["depmap_id","cell_line_display_name" ]):
                continue

            dataset_values = dataset[name].values
            dataset_values = pd.to_numeric(dataset_values, errors='coerce')

            # Check if the data is constant
            if len(np.unique(col_data)) == 1 or len(np.unique(dataset_values)) == 1:
                continue

            # Identify non-NaN values in both columns
            valid_mask = ~np.isnan(col_data) & ~np.isnan(dataset_values)
            valid_count = np.sum(valid_mask)  # Count the number of valid values

            # Calculate Spearman correlation and p-value
            spearman_corr_value, spearman_p_value = stats.spearmanr(
                dataset_values[valid_mask], col_data[valid_mask], nan_policy='omit'
            )

            # Append the results to the list
            results.append([name, col_name,valid_count, spearman_corr_value, spearman_p_value])

        # Print progress every 2 features
        if (list(features.columns).index(col_name) + 1) % 2 == 0 or col_name == features.columns[-1]:
            print(f"Number of features done: {list(features.columns).index(col_name) + 1}")

    # Convert the collected results to a DataFrame
    final_result_df = pd.DataFrame(
        results,
        columns=['Gene', 'Feature Name', 'Valid Count', 'Spearman Correlation', 'Spearman P-value']
    )

    # Extract p-values
    p_values = final_result_df['Spearman P-value']
    

    # Calculate FDR Adjusted P-values using Benjamini-Hochberg
    _, fdr_adjusted_p_values, _, _ = multipletests(p_values, method='fdr_bh')

    # Add the FDR Adjusted P-values to the DataFrame
    final_result_df['FDR Adjusted P-value'] = fdr_adjusted_p_values

    # Save the final DataFrame to a CSV file
    final_result_df.to_csv('gene_expression_correlations.csv', index=False)

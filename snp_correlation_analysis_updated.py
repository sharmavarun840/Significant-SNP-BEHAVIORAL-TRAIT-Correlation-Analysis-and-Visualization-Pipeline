
import pandas as pd
import scipy.stats as stats
import seaborn as sns
import matplotlib.pyplot as plt

def load_data(file_path, sheet_name='Sheet3'):
    # Load the Excel file
    xls = pd.ExcelFile(file_path)
    df = pd.read_excel(xls, sheet_name=sheet_name)
    return df

def preprocess_data(df):
    # Select relevant columns for analysis
    columns_of_interest = ['Physical Activity Level code', 'Outside food Code', 'Strees profile'] +                           [col for col in df.columns if col.startswith('rs')]
    # Filter and clean the dataset
    cleaned_df = df[columns_of_interest].dropna()
    return cleaned_df

def calculate_correlations(df, codes):
    # One-hot encode the genetic variant columns
    encoded_df = pd.get_dummies(df, columns=[col for col in df.columns if col.startswith('rs')])
    # Calculate correlation coefficients
    correlation_matrix = encoded_df.corr()
    correlations = correlation_matrix.loc[codes, [col for col in encoded_df.columns if col.startswith('rs')]]
    return correlations.T

def calculate_p_values(df, codes):
    p_values = pd.DataFrame(index=df.columns, columns=codes)
    for code in codes:
        for col in df.columns:
            if col.startswith('rs'):
                corr_coeff, p_value = stats.pearsonr(df[code], df[col])
                p_values.at[col, code] = p_value
    return p_values

def filter_significant_correlations(correlations, p_values, threshold=0.05):
    # Combine correlations and p-values, and filter for significance
    results_df = correlations.copy()
    results_df.columns = [f'Correlation with {col}' for col in results_df.columns]
    for col in correlations.columns:
        results_df[f'P-value for {col}'] = p_values[col]
    
    significant_correlations = results_df[
        (results_df['P-value for Physical Activity Level code'] < threshold) |
        (results_df['P-value for Outside food Code'] < threshold) |
        (results_df['P-value for Strees profile'] < threshold)
    ]
    return significant_correlations

def visualize_significant_correlations(significant_correlations):
    # Adjust SNP names and column labels
    significant_correlations.index = significant_correlations.index.str.split('_').str[0]
    significant_correlations = significant_correlations.rename(columns={
        'Correlation with Physical Activity Level code': 'Correlation with Physical Activity',
        'Correlation with Outside food Code': 'Correlation with Unhealthy Diet'
    })
    # Plot the heatmap
    plt.figure(figsize=(10, 6))
    sns.heatmap(significant_correlations.filter(like='Correlation'), annot=True, cmap='coolwarm', center=0)
    plt.title('Significant Correlations (P-value < 0.05)')
    plt.show()

def main():
    file_path = 'clean_data_cases_CoGSI.xlsx'  # Path to your file
    df = load_data(file_path)
    cleaned_df = preprocess_data(df)
    
    # Define the codes to correlate with
    behavioral_codes = ['Physical Activity Level code', 'Outside food Code', 'Strees profile']
    
    correlations = calculate_correlations(cleaned_df, behavioral_codes)
    p_values = calculate_p_values(cleaned_df, behavioral_codes)
    
    significant_correlations = filter_significant_correlations(correlations, p_values)
    visualize_significant_correlations(significant_correlations)

if __name__ == "__main__":
    main()

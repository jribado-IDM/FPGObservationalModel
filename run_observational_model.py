import pandas as pd
import numpy as np
import sys
import ast
import math
import os
import inspect
import json
from os.path import basename

from unified_sampling import run_sampling_model
from unified_metric_calculations import register_matrix, run_time_summaries, generate_het_barcode


#####################################################################################
# Utility functions
#####################################################################################
def get_default_config():
    """Return the default observational model configuration."""
    observational_model_config = {
        'hard_filters': {
            'symptomatics_only': False, 
            'monogenomic_infections_only': False,
            'day_snapshot': False
        },
        'intervention_start_month': 29, # Provide month where an intervention is applied. Currently any sampling pre/post intervention for a single intervention is supported. 
        'sampling_configs': {
            'random': {
                'method': 'random',
                'n_samples_year': 20,
                'replicates': 2,
                'method_params': {
                    'population_proportions': [1, 0], # Use to sample from the source or sink only, equally, etc. Within population comparisons of genetic metrics can be specified below - just make sure to total number of samples per year * proportion reflects the numbers you want per population.
                    'monogenomic_proportion': 0.2, # Set to False if sampling randomly 
                    'equal_monthly': False}
            },
            # 'seasonal': {
            #     'method': 'seasonal',
            #     'n_samples_year': 20,
            #     'replicates': 2,
            #     'method_params': {
            #         'season': 'full', # Options: full or peak; currently hardcoded to match Senegal's seasonality; update for other scenarios in unified_sampling.py
            #     }
            # },
            # 'age': { # Example of how to set-up a sampling scheme based on age, to mirror biased sampling such as school surveys and health facility comparisons. 
            #     'method': 'age',
            #     'n_samples_year': 15,
            #     'replicates': 1
            # }
            
        },
        'metrics': {
            'cotransmission_proportion': True,
            'complexity_of_infection': True, # Will calculate both true COI and effective COI from the unique strains identified in the sample.
            'heterozygosity': True,
            'identity_by_descent': False,
            'identity_by_state': True,
            'individual_ibx': True,
            'monogenomic_proportion': True,
            'rh': True,
            'unique_genome_proportion': True # Will calculate both the proportion of unique genomes in the sampled infections to replicate phasing and from monogenomic samples with an effective COI of 1  only to match barcode limits.
        },
        'subpopulation_comparisons': { # Supported for yearly and seasonal temporal sampling schemes, not age-based sampling. 
            'populations': True,  # Defined by the population node in EMOD
            'polygenomic': True,  # Is polygenomic = 1, else monogenomic = 0
            'symptomatic': True,  # Is symptomatic = 1, else asymptomatic = 0
            'age_bins': True      # Default age bins: 0-5, 5-15, 15+
        }
    }
    
    return observational_model_config


def update_matrix_indices(sample_df):
    """ 
    Update the sample_df DataFrame to include a new column 'original_nid' that stores the original maps the 'recursive_nid' values to a global order based on their unique values.
    """
    df = sample_df.copy()  # Fixed: use sample_df instead of undefined df
    
    # Parse the recursive_nid column
    df['original_nid'] = df['recursive_nid'].copy()
    df['original_nid'] = df['original_nid'].apply(ast.literal_eval)

    # Step 1: Get all unique recursive_nid values across all rows
    all_nids = []
    for nid_list in df['original_nid']:
        all_nids.extend(nid_list)

    # Get unique values and sort them
    unique_nids = sorted(set(all_nids))

    # Step 2: Create a mapping from nid to its global order
    nid_to_order = {nid: i for i, nid in enumerate(unique_nids)}

    # Step 3: Apply the mapping to create the order column
    def map_to_global_order(nid_list):
        return [nid_to_order[nid] for nid in nid_list]

    df['recursive_nid'] = df['original_nid'].apply(map_to_global_order)

    return df


def extract_sampled_infections(sample_df):
    """
    Remove rows where all sampling columns ('rep' columns) are NaN.
    
    Parameters:
        sampled_df (pd.DataFrame): DataFrame with sampling columns
        
    Returns:
        pd.DataFrame: DataFrame with rows removed where all sampling columns are NaN
    """
    # String to match in column names
    match_string = 'rep'
    
    # Identify columns matching the string
    sample_cols = sample_df.filter(regex=match_string)
    
    if sample_cols.empty:
        print("No 'rep' columns found - returning original DataFrame")
        return sample_df
    
    # Find rows where ALL sampling columns are NaN
    rows_all_nan = sample_cols.isna().all(axis=1)
    rows_to_drop = rows_all_nan.sum()
    
    print(f"Found {len(sample_cols.columns)} sampling columns: {sample_cols.columns.tolist()}")
    print(f"Dropping {rows_to_drop} rows where all sampling columns are NaN")
    
    # Keep rows where at least one sampling column has a value
    df_cleaned = sample_df[~rows_all_nan]
    
    return df_cleaned


# Save dictionaries of IBx distributions
def save_json(filename: str, json_data: object, indent: int = None, separators=None) -> None:
    with open(filename, 'w') as handle:
        json.dump(json_data, handle, cls=NumpyEncoder, indent=indent, separators=separators)
    return

# Add NumpyEncoder class for JSON serialization
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, np.integer):
            return int(obj)
        elif isinstance(obj, np.floating):
            return float(obj)
        elif isinstance(obj, np.bool_):
            return bool(obj)
        return super(NumpyEncoder, self).default(obj)

def make_json_serializable(obj):
    """Convert nested dictionary with tuple keys to JSON-serializable format"""
    if isinstance(obj, dict):
        return {str(k): make_json_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, tuple):
        return str(obj)
    else:
        return obj        


#####################################################################################
# Functions to read in all files
#####################################################################################
def run_observational_model(
    sim_name,
    emod_output_path,
    config_path,
    output_path,
    verbose=True):

    # Identify running parameters from config file or use defaults
    if os.path.isfile(config_path):
        with open(config_path, 'r') as file:
            config = json.load(file)
        if verbose:
            print(f"Loaded config from {config_path}")
    else:
        if verbose:
            print("Input config unavailable, using default config to specify running modes." )
        config = get_default_config()

    # Read in infection data
    infection_df_path = f'{emod_output_path}/infIndexRecursive-genomes-df.csv'
    if os.path.exists(infection_df_path):
        infection_df = pd.read_csv(infection_df_path)
        if verbose:
            print(f"Loaded data from {infection_df_path}: {len(infection_df)} records") 
    else:
        if verbose:
            # print("Input file unavailable, creating sample data for testing...")
            print(f"Error: {infection_df_path} not found. Loading test data.")
        infection_df = pd.read_csv('test_data/test_fpg_infections.csv')   


    # Run sampling model
    sample_df = run_sampling_model(
        input_df = infection_df,  
        config = config,  
        intervention_start_month = config['intervention_start_month']
    )
    sample_df = extract_sampled_infections(sample_df)
    sample_df = update_matrix_indices(sample_df)

    # Identify additional file for metric calculations - mainly IBx
    # Memory-map the genotype files (doesn't load to RAM) for access later
    user_specified_ibx = []
    ibd_matrix = None
    ibs_matrix = None

    if config['metrics']['identity_by_descent']:
        user_specified_ibx.append('ibd')  # Fixed: changed from 'ibx' to 'ibd'
        root_matrix_path = f'{emod_output_path}/roots.npy'
        if os.path.exists(root_matrix_path):
            ibd_matrix = np.load(root_matrix_path, mmap_mode='r')   
            register_matrix('ibd_matrix', ibd_matrix)
        else:
            print(f"Warning: {root_matrix_path} not found, IBD calculations will be skipped") 

    if config['metrics']['identity_by_state'] or config['metrics'].get('heterozygosity', True) or config['metrics']['rh']:
        user_specified_ibx.append('ibs')
        genotype_matrix_path = f'{emod_output_path}/variants.npy'
        if os.path.exists(genotype_matrix_path):
            ibs_matrix = np.load(genotype_matrix_path, mmap_mode='r')
        else:
            print(f"Error: {genotype_matrix_path} not found. Loading test data.") 
            ibs_matrix = np.load("test_data/test_variants.npy", mmap_mode='r')
        register_matrix('ibs_matrix', ibs_matrix)

    if config['metrics'].get('heterozygosity', True) and ibs_matrix is not None:
        # Generate barcode with Ns for heterozygosity calculations
        sample_df['barcode_with_Ns'] = sample_df.apply(lambda row: generate_het_barcode(ibs_matrix, row['recursive_nid']), axis=1)

        sample_df['heterozygosity'] = sample_df['barcode_with_Ns'].apply(
            lambda x: x.count('N')/len(x) if isinstance(x, list) and len(x) > 0 else 0
        )

    # Run metric calculations
    sample_df['group_year'] = sample_df['intervention_year'].copy() if 'intervention_year' in sample_df.columns else sample_df['simulation_year'].copy() 
    sample_df['group_month'] = sample_df['intervention_month'].copy() if 'intervention_month' in sample_df.columns else sample_df['continuous_month'].copy() 
    all_summaries, all_infection_ibx, all_ibx_dist_dict = run_time_summaries(
        sample_df, 
        subpop_config=config['subpopulation_comparisons'],
        user_ibx_categories=user_specified_ibx
    )

    # Save outputs
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    summary_output_filepath = f'{output_path}/{sim_name}_FPG_ModelSummaries.csv'
    all_summaries.to_csv(summary_output_filepath, index=False)
    sample_output_filepath = f'{output_path}/{sim_name}_FPG_SampledInfections.csv'
    # Merge in individual IBx results for sampled infections
    sample_df_merged = sample_df.merge(all_infection_ibx, on='infIndex', how='left')
    sample_df_merged.to_csv(sample_output_filepath, index=False)

    save_ibx_distributions = True
    if save_ibx_distributions:
        serializable_dict = make_json_serializable(all_ibx_dist_dict)
        for ibx_category, dist_dict in serializable_dict.items():
            dist_output_filepath = f'{output_path}/{sim_name}_{ibx_category}_distributions.json'
            save_json(dist_output_filepath, dist_dict)
            if verbose:
                print(f"Saved {ibx_category} distributions to {dist_output_filepath}")

    if verbose:
        print(f"Saved summary output to {summary_output_filepath}")
        print(f"Saved sample output to {sample_output_filepath}")

    return all_summaries, sample_df


if __name__ == "__main__":
    # Example usage
    sim_name = "test_simulation"
    emod_output_path = "./test_data"  # Adjust as needed
    config_path = "./config.json"     # Adjust as needed
    output_path = "output"          
    
    try:
        summaries, samples = run_observational_model(
            sim_name=sim_name,
            emod_output_path=emod_output_path,
            config_path=config_path,
            output_path=output_path
        )
        print("Model run completed successfully!")
    except Exception as e:
        print(f"Error running model: {e}")

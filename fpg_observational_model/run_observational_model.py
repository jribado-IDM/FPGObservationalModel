import pandas as pd
import numpy as np
import sys
import ast
import math
import os
import inspect
import json
from os.path import basename

from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from fpg_observational_model.unified_sampling import run_sampling_model
from fpg_observational_model.unified_metric_calculations import register_matrix, run_time_summaries, generate_het_barcode


#####################################################################################
# Utility functions
#####################################################################################
def get_default_config():
    """Return the default observational model configuration."""
    observational_model_config = {
        'hard_filters': {
            'symptomatics_only': True, 
            'monogenomic_infections_only': False,
            'day_snapshot': False
        },
        'intervention_start_month': 29, # Provide month where an intervention is applied. Currently any sampling pre/post intervention for a single intervention is supported. 
        'sampling_configs': {
            # 'random': {
            #     'method': 'random',
            #     'n_samples_year': 100,
            #     'replicates': 2,
            #     'method_params': {
            #         'population_proportions': [1, 0], # Use to sample from the source or sink only, equally, etc. Within population comparisons of genetic metrics can be specified below - just make sure to total number of samples per year * proportion reflects the numbers you want per population.
            #         'monogenomic_proportion': False, # Set to False if sampling randomly 
            #         'equal_monthly': False}
            # },
            # 'seasonal': {
            #     'method': 'seasonal',
            #     'n_samples_year': 100,
            #     'replicates': 1,
            #     'method_params': {
            #         'season': 'full', # Options: full or peak; currently hardcoded to match Senegal's seasonality; update for other scenarios in unified_sampling.py
            #     }
            # } 
        },
        'metrics': {
            'cotransmission_proportion': True,
            'complexity_of_infection': True, # Will calculate both true COI and effective COI from the unique strains identified in the sample.
            'heterozygosity': True,
            'identity_by_descent': False,
            'identity_by_state': True,
            'individual_ibx': True,
            'fws': True,
            'monogenomic_proportion': True,
            'rh': True,
            'unique_genome_proportion': True # Will calculate both the proportion of unique genomes in the sampled infections to replicate phasing and from monogenomic samples with an effective COI of 1  only to match barcode limits.
        },
        'subpopulation_comparisons': { # Supported for yearly and seasonal temporal sampling schemes, not age-based sampling. 
            'add_monthly': False,  # Whether to add monthly comparisons in addition to yearly comparisons for temporal sampling schemes
            'polygenomic': True,  # Is polygenomic = 1, else monogenomic = 0
            'symptomatic': False,  # Is symptomatic = 1, else asymptomatic = 0
            'age_bins': False     # Default age bins: 0-5, 5-15, 15+
        }
    }
    
    return observational_model_config


def load_matrix_safely(file_path, max_retries=3, use_local_copy=True):
    """
    Safely load numpy matrix with memory mapping using raw reconstruction.
    """
    import os
    import numpy as np
    import tempfile
    
    if not os.path.exists(file_path):
        print(f"File not found: {file_path}")
        return None
    
    print(f"Loading matrix with mmap from: {file_path}")
    
    # Method 1: Try standard mmap first
    try:
        array = np.load(file_path, mmap_mode='r')
        print(f"Successfully loaded with standard mmap, shape: {array.shape}")
        return array
    except Exception as e1:
        print(f"Standard mmap failed: {e1}")
    
    # Method 2: Raw reconstruction with manual mmap
    try:
        print(f"Method 2: Raw reconstruction with mmap")
        with open(file_path, 'rb') as f:
            # Read header info
            magic = f.read(6)
            if magic != b'\x93NUMPY':
                raise ValueError("Not a numpy file")
            
            major, minor = f.read(2)
            
            if major == 1:
                header_len = np.frombuffer(f.read(2), dtype=np.uint16)[0]
            else:
                header_len = np.frombuffer(f.read(4), dtype=np.uint32)[0]
            
            header = f.read(header_len).decode('latin1')
            
            # Clean up header
            import ast
            import re
            header_clean = re.sub(r'\s+', ' ', header.strip())
            header_dict = ast.literal_eval(header_clean)
            
            shape = header_dict['shape']
            dtype = header_dict['descr']
            fortran_order = header_dict['fortran_order']
            
            # Calculate data offset
            data_offset = f.tell()
            
        print(f"Parsed - Shape: {shape}, Dtype: {dtype}, Data offset: {data_offset}")
        
        # Create memory-mapped array directly from file
        array = np.memmap(file_path, dtype=dtype, mode='r', 
                         offset=data_offset, shape=shape,
                         order='F' if fortran_order else 'C')
        
        print(f"Successfully created memmap with shape: {array.shape}")
        return array
        
    except Exception as e2:
        print(f"Mmap reconstruction failed: {e2}")
        
        # Fallback: Load into memory then create temp mmap file
        try:
            print("Fallback: Creating temporary mmap file")
            # First load the data using Method 3 from before
            with open(file_path, 'rb') as f:
                # ... same header parsing as above ...
                f.seek(data_offset)
                data = np.frombuffer(f.read(), dtype=dtype)
                
                if fortran_order:
                    array = data.reshape(shape, order='F')
                else:
                    array = data.reshape(shape, order='C')
            
            # Create a temporary file for memory mapping
            temp_file = tempfile.NamedTemporaryFile(delete=False, suffix='.npy')
            np.save(temp_file.name, array)
            temp_file.close()
            
            # Load as mmap from temp file
            mmap_array = np.load(temp_file.name, mmap_mode='r')
            
            # Store temp filename for cleanup (you'd need to handle this)
            mmap_array._temp_file = temp_file.name
            
            print(f"Created temporary mmap file: {temp_file.name}")
            return mmap_array
            
        except Exception as e3:
            print(f"Temp mmap fallback failed: {e3}")
    
    print(f"All mmap methods failed for {file_path}")
    return None


def update_matrix_indices(sample_df):
    """ 
    Update the sample_df DataFrame to include a new column 'original_nid' that stores the original maps the 'recursive_nid' values to a global order based on their unique values.
    """
    df = sample_df.copy()  # Fixed: use sample_df instead of undefined df
    
    # Parse the recursive_nid column
    df['original_nid'] = df['recursive_nid'].copy()
    df['original_nid'] = df['original_nid'].apply(
    lambda x: ast.literal_eval(x) if isinstance(x, str) else x)

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
        config_path=None,
        config=None,
        output_path=None,
        verbose=True):
    """
    Run the observational model with either a config file or config dictionary.

    Parameters:
        sim_name: Name of the simulation
        emod_output_path: Path to EMOD output files
        config_path: Path to config JSON file (optional)
        config: Configuration dictionary (optional)
        output_path: Directory to save outputs
        verbose: Whether to print verbose output

    Note: If both config_path and config are provided, config_path takes precedence.
          If neither is provided, default config is used.
          Missing parameters in provided config will be filled with default values.
          Extra parameters not in default config will trigger a warning.
    """

    # Helper function to check for unknown keys
    def check_unknown_keys(default_dict, user_dict, path="config"):
        """Recursively check for keys in user_dict that don't exist in default_dict."""
        unknown_keys = []
        for key, value in user_dict.items():
            if key not in default_dict:
                unknown_keys.append(f"{path}.{key}")
            elif isinstance(value, dict) and isinstance(default_dict.get(key), dict):
                # Recursively check nested dictionaries
                nested_unknown = check_unknown_keys(default_dict[key], value, f"{path}.{key}")
                unknown_keys.extend(nested_unknown)
        return unknown_keys

    # Helper function to deep merge dictionaries
    def deep_merge(default_dict, override_dict):
        """Recursively merge override_dict into default_dict, preserving defaults for missing keys."""
        result = default_dict.copy()
        for key, value in override_dict.items():
            if key in result and isinstance(result[key], dict) and isinstance(value, dict):
                result[key] = deep_merge(result[key], value)
            else:
                result[key] = value
        return result

    # Start with default config
    default_config = get_default_config()

    # Add a warning if both config_path and config are provided
    if config_path and config is not None:
        print("Warning: Both config_path and config dictionary provided. "
              "config_path will take precedence.")

    # Determine which config to use and merge with defaults
    user_config = None
    config_source = None

    if config_path and os.path.isfile(config_path):
        try:
            with open(config_path, 'r') as file:
                user_config = json.load(file)
            config_source = config_path
            if verbose:
                print(f"Loaded config from {config_path}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Error parsing JSON from {config_path}: {e}")
        except Exception as e:
            raise ValueError(f"Error reading config file {config_path}: {e}")
    elif config is not None:
        user_config = config
        config_source = "provided dictionary"
        if verbose:
            print("Using provided config dictionary")

    # Check for unknown keys and warn user
    if user_config is not None:
        unknown_keys = check_unknown_keys(default_config, user_config)
        if unknown_keys:
            warning_msg = f"\nWARNING: Found unknown configuration parameters in {config_source}:"
            for key in unknown_keys:
                warning_msg += f"\n  - {key}"
            warning_msg += "\n\nThese parameters will be ignored. Please check for typos or refer to default config."
            warning_msg += "\nValid top-level parameters are: " + ", ".join(default_config.keys())
            print(warning_msg)

            # Optional: Ask user to confirm if they want to continue
            if verbose:
                print("\nContinuing with valid parameters merged with defaults...\n")

        # Merge user config with defaults
        config = deep_merge(default_config, user_config)
        if verbose:
            print("Merged user config with default values for missing parameters")
    else:
        config = default_config
        if verbose:
            print("No config provided, using default config to specify running modes.")

    # Read in infection data
    infection_df_path = f'{emod_output_path}/infIndexRecursive-genomes-df.csv'
    if os.path.exists(infection_df_path):
        infection_df = pd.read_csv(infection_df_path)
        if verbose:
            print(f"Loaded data from {infection_df_path}: {len(infection_df)} records")
    else:
        if verbose:
            print(f"Error: {infection_df_path} not found. Loading test data.")
        infection_df = pd.read_csv('test_data/test_fpg_infections.csv')

    # Run sampling model
    print(f"Config paramters for sampling model:\n {config}")    
    sample_df = run_sampling_model(
        input_df=infection_df,
        config=config,
        intervention_start_month=config['intervention_start_month']
    )
    sample_df = extract_sampled_infections(sample_df)
    sample_df['original_nid'] = sample_df['recursive_nid'].copy()
    
    # NOTE: Commented out to preserve original_nid without changing recursive_nid. In theory this updated index can be used to read in smaller genotype matrix below of infections sampled across all sampling schemes. Will require changes to functions that use the original_nid to map genomes to infections. 
    # sample_df = update_matrix_indices(sample_df)
    # matrix_indices = sample_df['recursive_nid'].tolist()
 
    # Identify additional file for metric calculations - mainly IBx
    # Memory-map the genotype files (doesn't load to RAM) for access later
    user_specified_ibx = []
    ibd_matrix = None
    ibs_matrix = None
    # Optional - included if need to filter out non-variant tracked sites, i.e. immunity markers or drugR for calculating genetic metrics only on neutral variant sites.
    variant_indices = None
    # variant_indices = []

    if config['metrics']['identity_by_descent']:
        user_specified_ibx.append('ibd')
        root_matrix_path = f'{emod_output_path}/roots.npy'
        if os.path.exists(root_matrix_path):
            ibd_matrix = load_matrix_safely(root_matrix_path)
            if variant_indices is not None and len(variant_indices) > 0:
                ibd_matrix = ibd_matrix[:, variant_indices]
            register_matrix('ibd_matrix', ibd_matrix)
        else:
            print(f"Warning: {root_matrix_path} not found, IBD calculations will be skipped")

    if config['metrics']['identity_by_state'] or config['metrics'].get('heterozygosity', True) or config['metrics'][
        'rh']:
        user_specified_ibx.append('ibs')
        genotype_matrix_path = f'{emod_output_path}/variants.npy'
     
        if os.path.exists(genotype_matrix_path):
            ibs_matrix = load_matrix_safely(genotype_matrix_path)
            if variant_indices is not None and len(variant_indices) > 0:
                ibs_matrix = ibs_matrix[:, variant_indices]
        else:
            print(f"Error: {genotype_matrix_path} not found. Loading test data.")
            ibs_matrix = np.load("../test_data/variants.npy", mmap_mode='r')
        register_matrix('ibs_matrix', ibs_matrix)

    if config['metrics'].get('heterozygosity', True) and ibs_matrix is not None:
        # Generate barcode with Ns for heterozygosity calculations
        sample_df['original_nid'] = sample_df['original_nid'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
        sample_df[['genotype_coi', 'barcode_with_Ns', 'heterozygosity']] = sample_df.apply(lambda row: generate_het_barcode(ibs_matrix, row['original_nid']), axis=1, result_type='expand')

    # Run metric calculations
    all_summaries, all_infection_ibx, all_ibx_dist_dict = run_time_summaries(
        sample_df,
        subpop_config=config['subpopulation_comparisons'],
        user_ibx_categories=user_specified_ibx
    )

    # Save outputs
    if output_path is None:
        output_path = "output"

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    summary_output_filepath = f'{output_path}/{sim_name}_FPG_ModelSummaries.csv'
    all_summaries.to_csv(summary_output_filepath, index=False)
    sample_output_filepath = f'{output_path}/{sim_name}_FPG_SampledInfections.csv'
    # Merge in individual IBx results for sampled infections
    if not all_infection_ibx.empty:
        # remove the monthly samples if in the df to have a smaller sampling only infection file
        if 'month_rep0' in sample_df.columns:
            sample_df = sample_df.drop(columns=['month_rep0'])
            sample_df = extract_sampled_infections(sample_df)
        sample_df = sample_df.merge(all_infection_ibx, on='infIndex', how='inner')
    sample_df.to_csv(sample_output_filepath, index=False)

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

    return


#####################################################################################
# Single run test
#####################################################################################
# Single file test
if __name__ == "__main__":
    import argparse

    # Set up argument parser
    parser = argparse.ArgumentParser(description='Run the FPG observational model')
    parser.add_argument('--sim_name', type=str, default='test_simulation',
                        help='Name of the simulation (default: test_simulation)')
    parser.add_argument('--emod_output_path', type=str, default='./test_data',
                        help='Path to EMOD output directory (default: ./test_data)')
    parser.add_argument('--config_path', type=str, default=None,
                        help='Path to configuration JSON file (optional)')
    parser.add_argument('--output_path', type=str, default='output',
                        help='Directory to save outputs (default: output)')
    parser.add_argument('--verbose', action='store_true',
                        help='Enable verbose output')
    # config dictionary option is omitted for CLI simplicity. Use config file instead.

    # Parse arguments
    args = parser.parse_args()

    # Example usage
    print(f"Running observational model with:")
    print(f"  Simulation name: {args.sim_name}")
    print(f"  EMOD output path: {args.emod_output_path}")
    print(f"  Config path: {args.config_path if args.config_path else 'Using default config'}")
    print(f"  Output path: {args.output_path}")
    print(f"  Verbose: {args.verbose}")

    try:
        run_observational_model(
            sim_name=args.sim_name,
            emod_output_path=args.emod_output_path,
            config_path=args.config_path,
            output_path=args.output_path,
            verbose=args.verbose
        )
        print("\nModel run completed successfully!")
    except Exception as e:
        import traceback

        print(f"\nError running model: {e}")
        if args.verbose:
            print("\nFull traceback:")
            print(traceback.format_exc())


            
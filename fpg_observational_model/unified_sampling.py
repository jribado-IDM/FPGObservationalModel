import itertools
import numpy as np
import pandas as pd
import os, sys
import math
import ast
from os.path import join, dirname, basename, exists
from pathlib import Path
from typing import List
from ast import literal_eval


def parse_list(s):
    """Convert string representation of list to actual Python list."""
    try:
        return ast.literal_eval(s)
    except Exception:
        return []

def adjust_time_columns(df, intervention_start_month=None):
    """
    Add continuous_month column and identify intervention month if specified.
    
    Parameters:
      df (pd.DataFrame): Input DataFrame with 'year' and 'month' columns
      intervention_start_month (int): Continuous month when intervention starts (optional)
    
    Returns:
      pd.DataFrame: DataFrame with continuous_month added
    """
    df = df.copy()
    
    # Just note the intervention month if specified (for asterisk marking)
    if intervention_start_month is not None and intervention_start_month > 0:
        print(f"Intervention starts at continuous month {intervention_start_month}")
        df['intervention_month'] = df['continuous_month'] - intervention_start_month
        df['intervention_year'] = df['intervention_month']  // 12
    
    return df


def calculate_infection_metrics(df):
    """
    Calculate metrics for each infection before any sampling decisions.
    
    For each infection, calculates:
    1. true_coi: number of items in genome_ids
    2. effective_coi: number of unique items in genome_ids  
    3. cotx: cotransmission status (NA if COI==1, True if single bite event, False if multiple)
    4. Optionally notes intervention month for season marking
    
    Parameters:
      df (pd.DataFrame): Input FPG data
      intervention_start_month (int): Optional continuous month when intervention starts
    
    Returns:
      pd.DataFrame: Data with added infection metrics
    """
    print(f"Calculating infection metrics for {len(df)} infections...")
    
    df = df.copy()
    
    # Step 2: Parse genome_ids and calculate COI
    df["original_nid"] = df["recursive_nid"].apply(parse_list)
    df["true_coi"] = df["original_nid"].apply(len)
    df["effective_coi"] = df["original_nid"].apply(lambda x: len(set(x)))
    # Keep only the single identifiable genome ids per infections
    # Consideration for future expansion, for effective COI > 1, keep density information for both to use in heterozygosity calculations
    df["recursive_nid"] = df["original_nid"].apply(lambda x: list(set(x)))
    
    # Step 3: Parse bite_ids for cotransmission
    df["bite_ids"] = df["bite_ids"].apply(parse_list)
    
    # Step 4: Calculate cotransmission (cotx)
    def calc_cotx(row):
        if row['effective_coi'] == 1:
            return None  # NA for monogenomic
        elif len(set(row['bite_ids'])) == 1:
            return True  # Single bite event = cotransmission
        else:
            return False  # Multiple bite events = superinfection
    
    df["cotx"] = df.apply(calc_cotx, axis=1)
    
    return df


def apply_emod_filters(infection_df, 
    fever_filter=False, 
    monogenomic_filter=False, 
    day_filter=False,
    other_filters=None):
    """
    Apply initial filters to the data before any sampling.
    
    Parameters:
      infection_df (pd.DataFrame): Input DataFrame
      fever_filter (bool): If True, keep only fever cases; if False, keep only non-fever cases
      monogenomic_filter (bool): If True, keep only monogenomic; if False, keep only polygenomic
      coi_filter (str): 'true_coi' or 'effective_coi' - which COI column to use
      other_filters (dict): Additional filters to apply {column: value}
    
    Returns:
      pd.DataFrame: Filtered DataFrame
    """
    df = infection_df.copy()
    original_size = len(df)

    # Apply timescale filter to have individuals included at most once per month, so filter out the daily infections that my not be changing 
    if 'continuous_month' not in df.columns:
        df = convert_month(df)
    df = filter_emod_infections(df, duplicate_seed=123)
    
    # Apply fever filter
    if fever_filter:
        if 'fever_status' not in df.columns:
            raise ValueError("fever_status column not found but fever_filter requested")
        
        fever_value = 1 if fever_filter else 0
        df = df[df['fever_status'] == fever_value]
        fever_remaining = len(df)
        print(f"    Fever filter: {original_size} -> {fever_remaining} samples "
              f"({'fever cases' if fever_filter else 'non-fever cases'} only)")
    
    # Apply COI filter
    if monogenomic_filter:
        coi_col = 'effective_coi' 
        if coi_col not in df.columns:
            raise ValueError(f"{coi_col} column not found but monogenomic_filter requested. Run calculate_infection_metrics() first.")
        
            df = df[df[coi_col] == 1]
            print(f"    COI filter: {len(df)} monogenomic infections kept")

    if day_filter:
        if df['day'].min() <= int(day_filter) <= df['day'].max():
            df = df[df['day'] == int(day_filter)]
        else:
            raise ValueError(f"User specified day {day_filter} is out of range of the simulation days({df['day'].min()} to {df['day'].max()}).")

    # Apply other filters
    if other_filters:
        for column, value in other_filters.items():
            if column not in df.columns:
                print(f"    Warning: Column '{column}' not found, skipping filter")
                continue
            
            before_filter = len(df)
            df = df[df[column] == value]
            after_filter = len(df)
            print(f"    Filter {column}={value}: {before_filter} -> {after_filter} samples")
    
    if len(df) == 0:
        print("    Warning: All data filtered out!")
    
    return df


#####################################################################################
# Utility functions - Sampling
#####################################################################################
def assign_season_group(row):
    """Assign wet/dry season groups based on month. """
    
    year = int(row['simulation_year'])
    month = int(row['month'])
    
    if month == 1:
        return f"Wet season: {year-1}-08 to {year}-01"
    elif month >= 8:
        return f"Wet season: {year}-08 to {year+1}-01"
    else:
        return f"Dry season: {year}-02 to {year}-07"


def assign_peak_group(row):
    """Assign peak/off-peak season groups. """
    year_col = 'simulation_year' if 'simulation_year' in row else 'year'
    
    year = int(row['simulation_year'])
    month = int(row['month'])
    
    if month >= 10 and month <= 12:
        return f"Peak wet: {year}-10 to {year}-12"
    elif 3 <= month <= 6:
        return f"Peak dry: {year}-03 to {year}-06"
    else:
        return f"Off-peak"


def convert_month(df):
    """ Convert simulation year and month to a continuous variable. """
    if 'month' in df.columns and 'simulation_year' in df.columns:
        df['continuous_month'] = (df['simulation_year'] * 12) + df['month']
    elif 'month' in df.columns and 'year' in df.columns:
        df['continuous_month'] = (df['year'] * 12) + df['month']
    return df


def create_robust_random_seed(base_seed, replicate, year=None, population=None, extra=None):
    """Create a robust random seed that avoids correlations."""
    seed_components = [base_seed, replicate * 1000]
    if year is not None:
        seed_components.append(int(year) * 100)
    if population is not None:
        seed_components.append(population * 10)
    if extra is not None:
        seed_components.append(hash(str(extra)) % 1000)
    
    combined_seed = sum(seed_components) % (2**31 - 1)
    return combined_seed


def handle_insufficient_samples(available, requested, group_info=""):
    """Handle cases where requested samples exceed available samples."""
    if available < requested:
        print(f"    Warning {group_info}: Requested {requested} but only {available} available. Using all.")
        return available
    return requested


def n_samples_by_pop(infection_df, n_samples_year, population_proportions=False):
    """ 
    Returns the number of samples per population based on total samples and optional population fractions.
    """
    n_pop = len(infection_df['population'].unique())
    
    if n_pop == 1:
        return [n_samples_year]

    print(f"    Dataset contains {n_pop} populations.")
    
    if not population_proportions or len(population_proportions) == 0:
        samples_per_pop = int(np.floor(n_samples_year / n_pop))
        pop_samples = [samples_per_pop] * n_pop
        print(f"    Equal distribution: {samples_per_pop} samples per population")
    else:
        if len(population_proportions) != n_pop:
            raise ValueError(f"Population fractions length ({len(population_proportions)}) doesn't match number of populations ({n_pop})")
        
        if abs(sum(population_proportions) - 1.0) > 0.01:
            raise ValueError("Population fractions must sum to 1.0")
        
        pop_samples = [int(n_samples_year * frac) for frac in population_proportions]
        print(f"    Weighted distribution: {pop_samples} samples per population (fractions: {population_proportions})")
    
    return pop_samples


#####################################################################################
# Helper functions for subset operations
#####################################################################################
def filter_emod_infections(infection_df, 
    duplicate_window='continuous_month', 
    duplicate_seed=123,
    is_test=False):
    """Filter EMOD infections to avoid duplicates."""
    times = np.unique(infection_df.day.values) if 'day' in infection_df.columns else []
    if len(times) > 0 and sorted(times) == list(range(min(times), max(times)+1)):
        print("    Warning: Population is output per day, individuals may be sampled multiple times.\n"
              "    Choosing one infection per individual at random for each time window.")

    if duplicate_window == 'continuous_month':    
        print(f"    Picking one infection per individual per month.")
        if is_test:
            infections_sub = infection_df.sort_values(['infIndex','continuous_month'],ascending=True).groupby(['IndividualID', 'continuous_month']).head(1)
        else:
            infections_sub = infection_df.groupby(['IndividualID', 'continuous_month']).sample(n=1, random_state=duplicate_seed)
    else:
        year_col = 'intervention_year' if 'intervention_year' in infection_df.columns else 'simulation_year'
        if year_col == 'intervention_year':
            print(f"    Picking one infection per individual per intervention adjusted year")
        else:
            print(f"    Picking one infection per individual per simulation year.")

        if is_test:
            infections_sub = infection_df.sort_values(['infIndex', year_col],ascending=True).groupby(['IndividualID', year_col]).head(1)
        else:
            infections_sub = infection_df.groupby(['IndividualID', year_col]).sample(n=1, random_state=duplicate_seed)

    return infections_sub


def validate_subset_inputs(infection_df, n_samples_year, replicates, method_name):
    """Validate common inputs for subset functions."""
    if infection_df is None or infection_df.empty:
        raise ValueError(f"{method_name}: Input DataFrame is empty or None")
    
    if n_samples_year <= 0:
        raise ValueError(f"{method_name}: n_samples_year must be positive, got {n_samples_year}")
    
    if replicates <= 0:
        raise ValueError(f"{method_name}: replicates must be positive, got {replicates}")
    
    required_cols = ['population', 'simulation_year', 'infIndex']
    missing_cols = [col for col in required_cols if col not in infection_df.columns]
    if missing_cols:
        raise ValueError(f"{method_name}: Missing required columns: {missing_cols}")



#####################################################################################
# Main subset functions (corrected)
#####################################################################################
def subset_randomly(infection_df, 
                   n_samples_year, 
                   replicates,
                   scheme='random',
                   monogenomic_proportion=False,
                   equal_monthly=False, 
                   population_proportions=False, 
                   base_seed=418):
    """
    Randomly subset samples and add columns to the original dataframe.
    
    Parameters:
      infection_df (pd.DataFrame): Input DataFrame
      n_samples_year (int): Number of samples per year per population
      replicates (int): Number of replicates
      scheme (str): Base name for columns (e.g., 'random')
      equal_monthly (bool): Whether to sample equally across months
      population_proportions (list): Optional population fractions (must sum to 1)
      base_seed (int): Base random seed
      
    Returns:
      pd.DataFrame: Original dataframe with new sampling columns added
    """
    # Validate inputs
    validate_subset_inputs(infection_df, n_samples_year, replicates, "subset_randomly")
    
    # Start with copy of original dataframe
    result_df = infection_df.copy()
    
    # Get sampling numbers per population
    n_per_pop = n_samples_by_pop(infection_df, n_samples_year, population_proportions)
    populations = sorted(infection_df['population'].unique())
    
    # Process each replicate
    for rep in range(replicates):
        # Create column name for this replicate
        sample_col_name = f"{scheme}_{n_samples_year}_rep{rep+1}"
        # Initialize column with NA
        result_df[sample_col_name] = pd.NA
        # Track sampled indices for this replicate
        sampled_indices = []
        
        for pop_id, pop_n_samp in zip(populations, n_per_pop):
            pop_df = infection_df[infection_df['population'] == pop_id]
            
            for year in sorted(infection_df['group_year'].unique()):
                year_df = pop_df[pop_df['group_year'] == year]
                
                if year_df.empty:
                    continue
                
                available_samples = len(year_df)
                target_samples = handle_insufficient_samples(
                    available_samples, pop_n_samp, 
                    f"Pop {pop_id}, Year {year}"
                )
                
                if target_samples == 0:
                    continue
                
                seed = create_robust_random_seed(base_seed, rep, year, pop_id)

                if monogenomic_proportion:
                    if not isinstance(monogenomic_proportion, float):
                        monogenomic_proportion = float(monogenomic_proportion)

                    if 0 < monogenomic_proportion < 1:
                        # print("Subsetting samples based on a targeted ",  monogenomic_proportion,  "monogenomic proportion",)
                        polygenomic_proportion = 1-monogenomic_proportion

                        mono_take = min(len(year_df[year_df['effective_coi'] == 1]), int(target_samples * monogenomic_proportion))
                        poly_take = min(len(year_df[year_df['effective_coi'] > 1]), int(target_samples * polygenomic_proportion))

                        year_sampled = pd.concat(
                            [year_df[year_df['effective_coi'] == 1].sample(mono_take, random_state=seed),
                            year_df[year_df['effective_coi'] > 1].sample(poly_take, random_state=seed)]
                        )
                        sampled_indices.extend(year_sampled.index.tolist())
                
                elif equal_monthly and 'continuous_month' in year_df.columns:
                    # Equal monthly sampling
                    unique_months = sorted(year_df['continuous_month'].unique())
                    samples_per_month = n_samples_year // 12
                    
                    for month in unique_months:
                        month_df = year_df[year_df['continuous_month'] == month]
                        if month_df.empty:
                            continue
                        
                        month_available = len(month_df)
                        month_take = min(month_available, samples_per_month)
                        if month_take > 0:
                            month_seed = seed + month  # Simple seed modification
                            month_sampled = month_df.sample(month_take, random_state=month_seed)
                            sampled_indices.extend(month_sampled.index.tolist())
                else:
                    # Simple random sampling
                    year_sampled = year_df.sample(target_samples, random_state=seed)
                    sampled_indices.extend(year_sampled.index.tolist())
        
        # Mark sampled infections with 1
        result_df.loc[sampled_indices, sample_col_name] = 1
    
    return result_df


def subset_by_seasons(infection_df, n_samples_year, replicates,
                      scheme="seasonal", season='full', base_seed=418):
    """
    Subset samples based on seasonal groupings.
    """
    result_df = infection_df.copy()
    
    # Apply seasonal groupings
    if season == "full":
        result_df['season_group'] = result_df.apply(assign_season_group, axis=1)
    elif season == "peak":
        result_df['season_group'] = result_df.apply(assign_peak_group, axis=1)
    else:
        raise ValueError("Only season_groups \"full\" or \"peak\" allowed as options.")

    validate_subset_inputs(result_df, n_samples_year, replicates, "subset_by_seasons")   
    
    # Process each replicate
    for rep in range(replicates):
        # Create column name for this replicate
        sample_col_name = f"{scheme}{season.capitalize()}_{n_samples_year}_rep{rep+1}"
        
        # Initialize column with NA
        result_df[sample_col_name] = pd.NA
        
        # Track sampled data for this replicate
        sampled_data = {}  # {index: season_label}
        
        season_groups = result_df['season_group'].unique()
        for season_name in season_groups:
            season_df = result_df[result_df['season_group'] == season_name]
            
            # Each season gets exactly n_samples_year samples
            available_samples = len(season_df)
            target_samples = handle_insufficient_samples(
                available_samples, n_samples_year, 
                f"Season {season_name}"
            )
            if season_df.empty or target_samples == 0:
                continue
            
            take = min(target_samples, available_samples)
            if take > 0:
                seed = create_robust_random_seed(base_seed, rep, extra=hash(season_name))
                sampled = season_df.sample(take, random_state=seed)
                
                # Store the season label for each sampled index
                for idx in sampled.index:
                    sampled_data[idx] = season_name
        
        # Set the season labels for sampled infections
        for idx, season_label in sampled_data.items():
            result_df.loc[idx, sample_col_name] = season_label

    return result_df


def subset_by_age(infection_df, n_samples_year, replicates,
                  scheme="ageBins",
                  age_bins=None, 
                  age_bin_weights=None,
                  base_seed=418):
    """
    Subset samples based on age bins and add columns to the original dataframe.
    
    Parameters:
      infection_df (pd.DataFrame): Input DataFrame
      n_samples_year (int): Number of samples per year (or per age bin if weights provided)
      replicates (int): Number of replicates
      scheme (str): Base name for columns (e.g., 'ageBins')
      age_bins (list): Age bin boundaries in days [0, bin1, bin2, ..., max]
      age_bin_weights (list): Optional weights for age bins (must sum to 1)
      base_seed (int): Base random seed
      
    Returns:
      pd.DataFrame: Original dataframe with new age bin sampling columns added
    """
    # Validate inputs
    validate_subset_inputs(infection_df, n_samples_year, replicates, "subset_by_age")
    
    # Start with copy of original dataframe
    result_df = infection_df.copy()
    
    # Set up age bins
    if age_bins is None:
        days_per_year = 365.25
        age_bins = [0, int(days_per_year * 5), int(days_per_year * 15), int(result_df['age_day'].max() + 1)]
        age_bin_labels = ['0-5yrs', '5-15yrs', '15+yrs']
    else:
        age_bin_labels = [f"age_{age_bins[i]//365}-{age_bins[i+1]//365}yrs" for i in range(len(age_bins)-1)]
    
    # Create age bin column
    result_df['age_bin'] = pd.cut(result_df['age_day'], bins=age_bins, labels=age_bin_labels, include_lowest=True)
    print(f"Created age bins: {result_df['age_bin'].value_counts().to_dict()}")
    
    # Calculate samples per age bin
    if age_bin_weights is not None:
        if len(age_bin_weights) != len(age_bin_labels):
            raise ValueError(f"Age bin weights length doesn't match number of age bins")
        if abs(sum(age_bin_weights) - 1.0) > 0.01:
            raise ValueError("Age bin weights must sum to 1.0")
        
        samples_per_age_bin = {label: int(n_samples_year * weight) 
                              for label, weight in zip(age_bin_labels, age_bin_weights)}
    else:
        # Equal distribution across age bins
        samples_per_age_bin = {label: n_samples_year // len(age_bin_labels) 
                              for label in age_bin_labels}
        # Distribute remainder
        remainder = n_samples_year % len(age_bin_labels)
        for i, label in enumerate(age_bin_labels[:remainder]):
            samples_per_age_bin[label] += 1
    
    # Process each replicate
    for rep in range(replicates):
        # Create column name for this replicate
        col_name = f"{scheme}_{n_samples_year}_rep{rep+1}"
        
        # Initialize column with NA
        result_df[col_name] = pd.NA
        
        # Track sampled indices for this replicate
        sampled_indices = []
        
        # Sample from each age bin
        for age_bin_label in age_bin_labels:
            age_df = result_df[result_df['age_bin'] == age_bin_label]
            
            # Get target samples for this age bin
            target_samples = samples_per_age_bin[age_bin_label]
            if age_df.empty or target_samples == 0:
                continue
            
            available = len(age_df)
            take = min(target_samples, available)
            
            if take > 0:
                seed = create_robust_random_seed(base_seed, rep, extra=hash(age_bin_label))
                sampled = age_df.sample(take, random_state=seed)
                
                # Add sampled indices to list
                sampled_indices.extend(sampled.index.tolist())
        
        # Mark sampled infections with 1
        result_df.loc[sampled_indices, col_name] = 1
    
    return result_df


#####################################################################################
# Main dispatcher function
#####################################################################################
def run_sampling_functions(infection_df, sampling_config, **kwargs):
    """
    Enhanced main function to run any subset sampling method.
    
    Parameters:
      infection_df (pd.DataFrame): Input DataFrame
      sampling_config (dict): Configuration dictionary with method, n_samples_year, replicates, method_params
      **kwargs: Additional arguments (e.g., base_seed)
    
    Returns:
      pd.DataFrame: DataFrame with sampling columns added
    """
    
    method = sampling_config['method']
    n_samples_year = sampling_config['n_samples_year']
    replicates = sampling_config['replicates']
    method_params = sampling_config.get('method_params', {})
    
    valid_methods = ['random', 'seasonal', 'age']
    if method not in valid_methods:
        raise ValueError(f"Unknown subset method: {method}. "
                         f"Available methods: {valid_methods}")
    
    print(f"\n=== Running {method} sampling ===")
    print(f"  Input: {len(infection_df)} records, {n_samples_year} samples/year, {replicates} replicates")
    
    base_seed = kwargs.get('base_seed', 418)
    
    try:
        if method == 'random':
            return subset_randomly(
                infection_df, n_samples_year, replicates,
                scheme='random',
                equal_monthly=method_params.get('equal_monthly', False),
                monogenomic_proportion=method_params.get('monogenomic_proportion'),
                population_proportions=method_params.get('population_proportions'),
                base_seed=base_seed
            )
            
        elif method == 'seasonal':
            return subset_by_seasons(
                infection_df, n_samples_year, replicates,
                scheme='seasonal',
                season=method_params.get('season', 'full'),
                base_seed=base_seed
            )
            
        elif method == 'age':
            return subset_by_age(
                infection_df, n_samples_year, replicates,
                scheme='ageBins',
                age_bins=method_params.get('age_bins'),
                age_bin_weights=method_params.get('age_bin_weights'),
                base_seed=base_seed
            )
            
    except Exception as e:
        print(f"  Error in {method} sampling: {str(e)}")
        raise



def process_config_filters(config):
    """Convert config filter settings to function parameters."""
    hard_filters = config.get('hard_filters', {})
    
    fever_filter = None
    monogenomic_filter = None
    
    # Convert string values to appropriate types
    if hard_filters.get('symptomatics_only') not in [False, None]:
        fever_filter = True
    
    if hard_filters.get('monogenomic_infections_only') not in [False, None]:
        monogenomic_filter = True
    
    return {
        'fever_filter': fever_filter,
        'monogenomic_filter': monogenomic_filter,
        'other_filters': None  # Could be extended
    }


def run_sampling_model(input_df, config, intervention_start_month=None, verbose=True):
    """
    Main function to run the observational model sampling pipeline.
    
    Parameters:
        input_df (pd.DataFrame): Input FPG DataFrame (required)
        config (dict): Configuration dictionary (required)
        intervention_start_month (int, optional): Continuous month when intervention starts
        verbose (bool): Whether to print progress messages
    
    Returns:
        pd.DataFrame: Final dataframe with all sampling columns added
    
    Raises:
        ValueError: If input_df or config is None
    """
    if verbose:
        print("=== Running Observational Model Sampling Pipeline ===")
    
    # Validate required inputs
    if input_df is None:
        raise ValueError("input_df is required and cannot be None")
    
    if config is None:
        raise ValueError("config is required and cannot be None")
    
    df = input_df.copy()
    
    # Add simulation_year column if it doesn't exist
    if 'simulation_year' not in df.columns and 'year' in df.columns:
        df['simulation_year'] = df['year'].copy()    
    
    
    try:
        # Step 1: Apply hard filters
        if verbose:
            print("\n=== Step 1: Apply hard filters ===")
        filter_params = process_config_filters(config)
        df_filtered = apply_emod_filters(df, **filter_params)

        intervention_start_month = config.get('intervention_start_month', intervention_start_month)
        if intervention_start_month is not None or intervention_start_month > 0:
            df_filtered = adjust_time_columns(df_filtered, intervention_start_month=intervention_start_month)

        df_filtered['group_year'] = df_filtered['intervention_year'].copy() if 'intervention_year' in df_filtered.columns else df['simulation_year'].copy()
        if config['subpopulation_comparisons'].get('add_monthly'):
            df_filtered['group_month'] = df_filtered['intervention_month'].copy() if 'intervention_month' in df_filtered.columns else df['continuous_month'].copy() 

        # Step 2: Calculate infection metrics
        if verbose:
            print("\n=== Step 2: Calculate individual infection metrics ===")
        df_metrics = calculate_infection_metrics(df_filtered.copy())
        
        # Step 3: Run sampling for each method
        if verbose:
            print("\n=== Step 3: Run sampling methods ===")
        
        # Start with the base dataframe
        final_df = df_metrics.copy()


        # Apply each sampling method
        for sampling_name, sampling_config in config['sampling_configs'].items():
            if verbose:
                print(f"\n--- Processing {sampling_name} sampling ---")
            
            # Run the sampling method
            sampled_df = run_sampling_functions(df_metrics, sampling_config)

            if config['subpopulation_comparisons'].get('add_monthly'):
                sampled_df['month_rep0'] = 1 
            
            # Merge the sampling columns back to the main dataframe
            sampling_columns = [col for col in sampled_df.columns 
                              if sampling_name in col and 'rep' in col]
            
            for col in sampling_columns:
                final_df[col] = sampled_df[col]
        
        if config['subpopulation_comparisons'].get('add_monthly'):
            final_df['month_rep0'] = 1  # All infections included for monthly analysis
            print("Added 'month_rep0' column for monthly subpopulation comparisons")
        

        # Final summary
        if verbose:
            print(f"\n=== Final Results ===")
            print(f"Final dataframe shape: {final_df.shape}")
            
            # Show sampling column summary
            sampling_cols = [col for col in final_df.columns 
                            if any(method in col for method in ['random', 'seasonal', 'age', 'month']) and 'rep' in col]
            print(f"Sampling columns created: {sampling_cols}")
            
            for col in sampling_cols:
                sampled_count = final_df[col].notna().sum()
                print(f"  {col}: {sampled_count} samples")
        
        return final_df
        
    except Exception as e:
        print(f"Error in observational model pipeline: {str(e)}")
        raise

    
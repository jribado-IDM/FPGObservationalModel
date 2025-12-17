import pandas as pd
import numpy as np
import ast
from itertools import combinations, chain
from collections import Counter
from scipy import stats
from sklearn.linear_model import LinearRegression
import idm, tskit

#####################################################################################
# Helper functions to access variables from calling module
#####################################################################################
# Matrix registry to store matrices from calling module
_matrix_registry = {}

def register_matrix(name, matrix):
    """Register a matrix for IBx calculations"""
    _matrix_registry[name] = matrix
    print(f"Registered matrix: {name}")

def get_matrix(name):
    """Get a registered matrix"""
    if name in _matrix_registry:
        return _matrix_registry[name]
    else:
        available = list(_matrix_registry.keys())
        raise KeyError(f"Matrix '{name}' not found. Available: {available}")


#####################################################################################
# Group identification for calculations
#####################################################################################
def identify_nested_comparisons(df, sampling_column_name, 
    config = None, 
    add_monthly=False):
    """
    Generate a list of infections within sampling schemes for looping through nested comparisons. 
    """
    nested_indices = {}

    # Specifying time groups
    if 'seasonal' not in sampling_column_name:
        time_group = "group_year"
        nested_indices[time_group] = df.groupby(time_group)['infIndex'].apply(list).to_dict()

        if add_monthly:
            nested_indices['group_month'] = df.groupby([time_group, 'group_month'])['infIndex'].apply(list).to_dict()         

    if 'seasonal' in sampling_column_name: 
        time_group = sampling_column_name
        if len(df[sampling_column_name].unique()) > 1:
            nested_indices['season_bins'] = df.groupby(sampling_column_name)['infIndex'].apply(list).to_dict()
        else:
            print("User specified comparisons by season, but only one season found.") 

    if 'age' in sampling_column_name: 
        if len(df[sampling_column_name].unique()) > 1:
            nested_indices['age_bins'] = df.groupby(['group_year', sampling_column_name])['infIndex'].apply(list).to_dict()
        else:
            print("User specified nested comparisons by age bin, but only one age bin available in the sample subset.")               

    # Specifying non-time groups; i.e. subgroups 
    if 'age' not in sampling_column_name and config is not None:
        if config.get('populations', False):
            if len(df['population'].unique()) > 1:  # FIXED: added len()
                nested_indices['populations'] = df.groupby([time_group, 'population'])['infIndex'].apply(list).to_dict()
            else:
                print("User specified nested comparisons by population, but only one population is available.")

        if config.get('polygenomic', False):
            df = df.copy()
            df['is_polygenomic'] = df['effective_coi'].apply(lambda x: True if x > 1 else False) 
            polygenomic_vals = df['is_polygenomic'].unique()
            if True in polygenomic_vals and False in polygenomic_vals:
                nested_indices['polygenomic'] = df.groupby([time_group, 'is_polygenomic'])['infIndex'].apply(list).to_dict()
            else:
                print("User specified nested comparisons by monogenomic or polygenomic infections, but only one group available.")   

        if config.get('symptomatic', False):  
            if len(df['fever_status'].unique()) > 1:
                nested_indices['symptomatic'] = df.groupby([time_group, 'fever_status'])['infIndex'].apply(list).to_dict()  
            else:
                print("User specified nested comparisons by fever status, but only one fever status is available.")       

        if config.get('age_bins', False):
            days_per_year = 365.25
            age_bins = [0, int(days_per_year * 5), int(days_per_year * 15), int(df['age_day'].max() + 1)]
            age_bin_labels = ['0-5yrs', '5-15yrs', '15+yrs']     
            df['age_bin'] = pd.cut(df['age_day'], bins=age_bins, labels=age_bin_labels, include_lowest=True)
            
            if len(df['age_bin'].unique()) > 1:
                nested_indices['age_bins'] = df.groupby([time_group, 'age_bin'])['infIndex'].apply(list).to_dict()  
            else:
                available_age_group = df['age_bin'].unique()[0]
                print(f"User specified nested comparisons by age bins, but only one age group, {available_age_group}, is available.")        

    return nested_indices

#####################################################################################
# Summary statistics calculations
#####################################################################################
def comprehensive_group_summary(group):
    """
    Calculate comprehensive summary statistics for infection data.
    
    Returns mean, median, std, min, max for continuous variables.
    """
    if len(group) == 0:
        return _empty_comprehensive_summary()
    
    # Monogenomics and polygenomics counts and COI stats
    true_coi_stats = _comprehensive_stats(group['true_coi'], 'true_coi')
    effective_coi_stats = _comprehensive_stats(group['effective_coi'], 'effective_coi')

    # Genome ID analysis
    all_genome_stats = _analyze_genome_ids(group['recursive_nid'], "all_genomes")
    mono_genome_stats = _analyze_genome_ids(group[group['effective_coi'] == 1]['recursive_nid'], "mono_genomes")
    
    # Cotransmission and superinfection analysis
    poly_series =  group[group['true_coi'] > 1]['cotx']
    cotxn_counts = _analyze_binary_in_subset(poly_series, 'cotransmission')

    # Combine all stats
    result = pd.Series({
        'n_infections': len(group),
        **true_coi_stats,
        **effective_coi_stats,
        **all_genome_stats,
        **mono_genome_stats,
        **cotxn_counts
    })
    
    if 'genotype_coi' in group.columns:
        genotype_coi_stats = _comprehensive_stats(group['genotype_coi'], 'genotype_coi')

        result = pd.concat([result, genotype_coi_stats])

    return result


def _comprehensive_stats(series, prefix):
    """Calculate comprehensive statistics with given prefix."""
    poly_mask = series > 1
    count = poly_mask.sum()
    prop  = round(count / len(series), 3)
    
    stats = series.describe()
    return pd.Series({
        f'{prefix}_poly_count': count,
        f'{prefix}_poly_prop': round(prop, 3),
        f'{prefix}_mean': round(stats['mean'], 3),
        f'{prefix}_median': round(stats['50%'], 3),
        f'{prefix}_std': round(stats['std'], 3),
        f'{prefix}_min': round(stats['min'], 3),
        f'{prefix}_max': round(stats['max'], 3),
        f'{prefix}_q25': round(stats['25%'], 3),
        f'{prefix}_q75': round(stats['75%'], 3)
    })

def _analyze_binary_in_subset(series, prefix):
    """Analyze binary column within a subset defined by mask."""
    if len(series) == 0:
        return {f'{prefix}_count': 0, f'{prefix}_prop': np.nan}
    
    mask = series == 1 
    count = mask.sum()
    prop = count / len(series)
    
    return {f'{prefix}_count': count, f'{prefix}_prop': prop}    


def _analyze_genome_ids(series, prefix):
    """Analyze genome IDs to get total count, unique count, and proportion."""
    all_genome_ids = []
    for sublist in series:
        if isinstance(sublist, list):
            all_genome_ids.extend(sublist)
    
    total_genomes  = len(all_genome_ids)
    unique_genomes = len(set(all_genome_ids))
    unique_prop = round(unique_genomes / total_genomes, 3) if total_genomes > 0 else np.nan
    
    return {
        f'{prefix}_total': total_genomes,
        f'{prefix}_unique': unique_genomes,
        f'{prefix}_unique_prop': unique_prop
        }

def _empty_comprehensive_summary():
    """Return comprehensive summary with NaN/0 values for empty groups."""
    base_stats = {
        'n_infections': 0,
        'cotransmission_count': 0,
        'cotransmission_prop': np.nan,
    }
    
    # Add comprehensive stats for both COI measures
    for prefix in ['effective_coi', 'true_coi', 'genotype_coi']:
        for stat in ['poly_count', 'poly_prop', 'mean', 'median', 'std', 'min', 'max', 'q25', 'q75']:
            base_stats[f'{prefix}_{stat}'] = np.nan

        for prefix in ['all_genomes', 'mono_genomes']:
            for stat in ['total', 'unique', 'unique_prop']:
                base_stats[f'{prefix}_{stat}'] = 0 if 'count' in stat else np.nan

    
    return pd.Series(base_stats)


################################################################################
# IBx calculation functions - i.e. organizes the pairs 
################################################################################
def update_ibx_index(filter_df):
    """ 
    For year specific IBX calculations, update the recursive_nid to a global order based on their unique values.
    """
    # Step 1: Get all unique recursive_nid values across all rows
    filter_df = filter_df.copy()
    all_nids = []
    for nid_list in filter_df['recursive_nid']:
        all_nids.extend(nid_list)

    # Get unique values and sort them
    unique_nids = sorted(set(all_nids))

    # Step 2: Create a mapping from nid to its global order
    nid_to_order = {nid: i for i, nid in enumerate(unique_nids)}

    # Step 3: Apply the mapping to create the order column
    def map_to_global_order(nid_list):
        return [nid_to_order[nid] for nid in nid_list]

    filter_df['ibx_nid'] = filter_df['recursive_nid'].apply(map_to_global_order)

    return filter_df


def calculate_ibx_matrix(df, genotypes, intervals=None):
    
    # update indices for the hash table
    df = update_ibx_index(df)
    df = df.explode(['recursive_nid', 'ibx_nid'])

    genotypes = idm.align_data(genotypes)
    if intervals is not None:
        hash_ibx, hashes, mapping = idm.calculate_ibx(genotypes, intervals=intervals)
    else:
        hash_ibx, hashes, mapping = idm.calculate_ibx(genotypes)

    if len(df.ibx_nid) == genotypes.shape[0]:
        ibx_index = pd.DataFrame(zip(df.ibx_nid.tolist(), [mapping[i] for i in hashes]),
            columns =['ibx_nid', "ibx_index"])
    else:
        ibx_index = pd.DataFrame(list(enumerate([mapping[i] for i in hashes])),
            columns =['ibx_nid', "ibx_index"])
 
    df = pd.merge(df, ibx_index, on='ibx_nid').reset_index(drop=True)
    cols = ['infIndex', 'recursive_nid', 'ibx_nid', 'ibx_index']
    hash_df = df[cols].drop_duplicates().reset_index(drop=True)
    
    return hash_df, hash_ibx


def ibx_distribution(indices, hash_ibx):
    """
    Fxn returns counts of pairwise values per group.
    Key: IBx value, value: counts
    """
    max_scaler = np.max(hash_ibx)
    ibx_dict = {}
    pairwise_hash = combinations(indices, 2)
    pairwise_counts = Counter(pairwise_hash) 
    # loop through pairwise combinations
    for key in pairwise_counts:
        weight = pairwise_counts[key]
        if len(set(key)) == 1:
            ibs = 1
        else:
            a, b = key  
            ibs = hash_ibx[a, b]/max_scaler
        # update the dictionary object
        ibs = round(ibs, 2)
        if ibs in ibx_dict:
            new_value = weight + ibx_dict[ibs]
            ibx_dict.update({ibs: new_value})
        else:
            ibx_dict.update({ibs: weight})
    return ibx_dict


################################################################################
# IBx summary functions - i.e. calculates summary statistics from the counts distribution dictionary 
################################################################################
# Matching the .describe() function with a dictionary input to update weighted mean and median functions used in the original version of the observational model
# Adapted with help with Claude Sonnet 4.0
def weighted_describe_scipy(summary_dict, ibx_prefix):
    """Calculate stats by expanding the weighted dictionary"""
    if not summary_dict:
        return pd.DataFrame()
    
    # Expand the dictionary to a list
    expanded_values = []
    for value, count in summary_dict.items():
        expanded_values.extend([value] * int(count))
    
    expanded_values = np.array(expanded_values)
    
    # Now use standard numpy/pandas functions
    count = len(expanded_values)
    mean = np.mean(expanded_values)
    std = np.std(expanded_values, ddof=1)  # Use ddof=1 for sample std to match pandas
    
    summary_data = {
        f'{ibx_prefix}_count': int(count),
        f'{ibx_prefix}_mean': round(mean, 3),
        f'{ibx_prefix}_std': round(std, 3),
        f'{ibx_prefix}_min': round(np.min(expanded_values), 3),
        f'{ibx_prefix}_25%': round(np.percentile(expanded_values, 25), 3),
        f'{ibx_prefix}_50%': round(np.median(expanded_values), 3),
        f'{ibx_prefix}_75%': round(np.percentile(expanded_values, 75), 3),
        f'{ibx_prefix}_max': round(np.max(expanded_values), 3)
    }
    
    return pd.DataFrame([summary_data]) 


#####################################################################################
# Run summaries and IBx calculations for nested groups
#####################################################################################
def process_nested_summaries(nested_indices, sampling_df, comprehensive_group_summary):
   
    summary_stats_list = []
    
    def add_summary(indices, comparison_type, year_group, subgroup=None):
        group_subset = sampling_df[sampling_df['infIndex'].isin(indices)]
        summary = comprehensive_group_summary(group_subset)
        
        if isinstance(summary, pd.Series):
            summary_dict = summary.to_dict()
        else:
            summary_dict = summary
            
        summary_dict.update({
            'comparison_type': comparison_type,
            'year_group': str(year_group),
            'subgroup': str(subgroup) if subgroup is not None else None
        })
        summary_stats_list.append(summary_dict)
    
    for comparison_type, data in nested_indices.items():
        if comparison_type in ['group_year', 'group_month']:
            for key, indices in data.items():
                if isinstance(key, tuple):
                    year_group = f"{key[0]}_{key[1]}" if comparison_type == 'group_month' else str(key[0])
                    add_summary(indices, comparison_type, year_group)
                else:
                    add_summary(indices, comparison_type, str(key))
        elif comparison_type.startswith('seasonal'):
            for key, indices in data.items():
                add_summary(indices, comparison_type, str(key))
        else:
            for key, indices in data.items():
                if isinstance(key, tuple) and len(key) == 2:
                    year_group, subgroup = key
                    if not isinstance(indices, float):  # Skip NaN values
                        add_summary(indices, comparison_type, str(year_group), str(subgroup))
                elif isinstance(indices, dict):
                    year_group = str(key)
                    for subgroup, sub_indices in indices.items():
                        if not isinstance(sub_indices, float):
                            add_summary(sub_indices, comparison_type, year_group, str(subgroup))
                else:
                    if not isinstance(indices, float):
                        add_summary(indices, comparison_type, str(key), None)
    
    if summary_stats_list:
        result_df = pd.DataFrame(summary_stats_list)
        print(f"Created summary DataFrame with shape: {result_df.shape}")
        return result_df
    else:
        print("Warning: No summary statistics generated")
        return pd.DataFrame()


def inf_ibx_summary(ibx_matrix, ibx_indices):
    """
    Run the IBx summary for a list of genome indices for each polygenomic infection.
    """
    distribution = ibx_distribution(ibx_indices, ibx_matrix)
    summary_stats = weighted_describe_scipy(distribution, "ibx")

    return summary_stats.iloc[0].to_dict()
    

def process_nested_ibx(df, gt_matrix, nested_indices, 
ibx_prefix,
individual_ibx_calculation=True,
save_ibx_distributions=True):
    """
    Calculate IBx for nested comparison groups.
    
    Args:
        df: DataFrame linking infection information to genotype indices
        gt_matrix: Matrix of genotypes, roots or alleles. 
        nested_indices: Dictionary with comparison types as keys and nested group data as values
        save_ibx_distributions: Option to save the dictionary of value with counts for the full pairwise distribution
        individual_ibx_calculation: Whether to calculate individual IBx values
        
    Returns:
        DataFrame with summary statistics for each group/subgroup
    """

    if 'group_year' in nested_indices.keys():
        all_year_indices = nested_indices['group_year']

    if 'season_bins' in nested_indices.keys():
        all_year_indices = nested_indices['season_bins']

    ibx_dist_dict, individual_ibx_dict = {}, {}
    ibx_summ_list = []
    for year_key, indices in all_year_indices.items():
        # Handle both string keys and tuple keys for year
        year = str(year_key) if not isinstance(year_key, tuple) else str(year_key[0])
        
        year_subset = df[df['infIndex'].isin(indices)]
        year_subset = update_ibx_index(year_subset) 

        # Step 1: Run pairwise IBx calculations once per year
        genome_indices = []
        for idx_list in year_subset['recursive_nid']:
            if isinstance(idx_list, list):
                genome_indices.extend(idx_list)

        matrix = get_matrix(gt_matrix)[genome_indices, :]
        print("Genotype matrix shape:", matrix.shape)
        ibx_indices, ibx_matrix = calculate_ibx_matrix(year_subset, matrix)

        # Add column with the ibx_index for each infection
        ibx_mapping = dict(zip(ibx_indices['ibx_nid'], ibx_indices['ibx_index']))
        year_subset['ibx_index'] = year_subset['ibx_nid'].apply(lambda nid_list: [ibx_mapping[nid] for nid in nid_list] if isinstance(nid_list, list) else None)
        
        # Step 2: Run IBx summaries for nested groups within each year
        for comparison_group, group_data in nested_indices.items():
            if comparison_group in ['group_year', 'season_bins']:
                # Simple year-level calculation
                indices = list(chain.from_iterable(year_subset['ibx_index'].tolist())) 

                if isinstance(indices, list) and len(indices) > 1:
                    distribution = ibx_distribution(indices, ibx_matrix)
                    summary_stats = weighted_describe_scipy(distribution, ibx_prefix)
                    
                    # Add metadata columns
                    result_row = summary_stats.iloc[0].to_dict()
                    result_row['comparison_type'] = comparison_group
                    result_row['year_group'] = year
                    result_row['subgroup'] = None
                    ibx_summ_list.append(result_row)

                    if save_ibx_distributions:
                        if comparison_group not in ibx_dist_dict:
                            ibx_dist_dict[comparison_group] = {}
                        # FIXED: Remove the if condition, just assign
                        ibx_dist_dict[comparison_group][year] = distribution

            else:
                # Handle nested groups
                for key, nested_data in group_data.items():
                    if isinstance(key, tuple):
                        key_year, subgroup = key
                        if str(key_year) == year:
                            if isinstance(nested_data, list) and len(nested_data) > 1:
                                subset_df = year_subset[year_subset['infIndex'].isin(nested_data)]
                                if not subset_df.empty:
                                    subgroup_ibx_indices = list(chain.from_iterable(subset_df['ibx_index'].tolist()))

                                    if subgroup_ibx_indices and len(subgroup_ibx_indices) > 1:
                                        distribution = ibx_distribution(subgroup_ibx_indices, ibx_matrix)
                                        summary_stats = weighted_describe_scipy(distribution, ibx_prefix)
                                        
                                        result_row = summary_stats.iloc[0].to_dict()
                                        result_row['comparison_type'] = comparison_group
                                        result_row['year_group'] = year
                                        result_row['subgroup'] = str(subgroup)
                                        ibx_summ_list.append(result_row)

                                        if save_ibx_distributions:
                                            if comparison_group not in ibx_dist_dict:
                                                ibx_dist_dict[comparison_group] = {}
                                            ibx_dist_dict[comparison_group][key] = distribution
            
        # Step 3: Individual IBx calculations for polygenomic infections
            if individual_ibx_calculation:
                polygenomic_subset = year_subset[year_subset['effective_coi'] > 1]
                if not polygenomic_subset.empty:
                    polygenomic_dict = dict(zip(polygenomic_subset['infIndex'], polygenomic_subset['ibx_index']))

                    for inf_id, ibx_list in polygenomic_dict.items():
                        distribution = ibx_distribution(ibx_list, ibx_matrix)
                        individual_ibx_dict[inf_id] = weighted_describe_scipy(distribution, ibx_prefix) 

    if ibx_summ_list:
        ibx_results_df = pd.DataFrame(ibx_summ_list)
    else:
        ibx_results_df = pd.DataFrame()

    if individual_ibx_dict:  
        individual_ibx_df = pd.concat(individual_ibx_dict, names=['infIndex', 'row_id']).reset_index(level=0)
    else:
        individual_ibx_df = pd.DataFrame()

    return ibx_results_df, individual_ibx_df, ibx_dist_dict
               

def run_time_summaries(sample_df,
subpop_config = None,
add_monthly = False, 
user_ibx_categories = None,
individual_ibx_calculation=True,
fws_calculation=True,
rh_calculation=True,
save_ibx_distributions=True):
    
    df = sample_df
    sampling_columns = df.filter(regex = "rep").columns.to_list()
    
    all_summary_dataframes = []
    all_rh_dataframes = []
    all_inf_ibx, all_inf_rh = [], []
    all_ibx_dist_dict = {}

    for sampling_column in sampling_columns:
        sampling_df = df[df[sampling_column].notna()]

        nested_dict = identify_nested_comparisons(sampling_df, sampling_column, config = subpop_config, add_monthly=add_monthly)
        print(f"Nested comparisons for {sampling_column}: {list(nested_dict.keys())}")

        # Get base summary statistics
        summary_stats = process_nested_summaries(nested_dict, sampling_df, comprehensive_group_summary)
        summary_stats.insert(0, 'sampling_scheme', sampling_column)

        # Process Fws calculations if specified
        if fws_calculation:
            fws_stats = process_nested_fws(nested_dict, sampling_df)
            summary_stats = summary_stats.merge(fws_stats, on=['comparison_type', 'year_group', 'subgroup'], how='left')

        # Process IBx categories if they exist
        if user_ibx_categories and len(user_ibx_categories) > 0: 
            for ibx_category in user_ibx_categories:
                ibx_dist_dict = {}
                print(f"\nProcessing IBx category: {ibx_category}")
                try:
                    ibx_summary, ibx_inf, ibx_dist_dict = process_nested_ibx(
                        sampling_df,  
                        f'{ibx_category}_matrix', 
                        nested_dict, 
                        ibx_prefix=ibx_category,
                        individual_ibx_calculation=individual_ibx_calculation,
                        save_ibx_distributions=save_ibx_distributions
                    )

                    if rh_calculation and ibx_category == 'ibs':
                        if 'polygenomic' not in nested_dict.keys():
                            if len(sampling_df['effective_coi'].unique()) != 1:
                                print("Warning: No polygenomic subpopulation comparisons found. Rerun IBx with polygenomic as a subpopulation comparisons.")
                            else:
                                continue    
                        
                        else:
                            # Only run R_h for year or seasonal level comparisons - can expand to do this with all subpopulations later if needed.

                            # Get all the distribution dictionaries for H_mono bootstrap calculations
                            monogenomic_dict = {k[0]:v for k,v in ibx_dist_dict['polygenomic'].items() if k[1]==False}

                            subpopulation_keys = [x for x in nested_dict.keys() if x in ['group_year', 'season_bins']]

                            for key in subpopulation_keys:
                                subpopulation_dict = nested_dict[key] 
                                yearly_rh_df = pd.DataFrame()
                                for sub_key, sub_data in subpopulation_dict.items():
                                    year_df = sampling_df[sampling_df['infIndex'].isin(sub_data)]  
                                    year_df = year_df.merge(pd.DataFrame(ibx_inf), on='infIndex', how='left')

                                    if sub_key in monogenomic_dict.keys():
                                        rh_summary, sample_rh = calculate_population_rh(year_df, monogenomic_dict[sub_key])

                                        rh_summary['comparison_type'] = key  
                                        rh_summary['year_group'] = str(sub_key)
                                        rh_summary['subgroup'] = None
                                        yearly_rh_df = pd.concat([yearly_rh_df, rh_summary], ignore_index=True)

                                        # Rename columns to indicate the sampling scheme for individual infections; will be used to confirm partner Rh groupings
                                        new_columns = [
                                            f"{sampling_column}-{col}" if 'rh' in col else col
                                            for col in sample_rh.columns
                                        ]
                                        sample_rh.columns = new_columns
                                        all_inf_rh.append(sample_rh)

                                summary_stats = summary_stats.merge(yearly_rh_df, on=['comparison_type', 'year_group', 'subgroup'], how='left')

                    if ibx_dist_dict:
                        if ibx_category not in all_ibx_dist_dict:
                            all_ibx_dist_dict[ibx_category] = {}
                        all_ibx_dist_dict[ibx_category][sampling_column] = ibx_dist_dict

                    if not ibx_inf.empty:
                        all_inf_ibx.append(ibx_inf)

                    if not ibx_summary.empty:
                        merge_keys = ['comparison_type', 'year_group', 'subgroup']
                        available_keys = [key for key in merge_keys if key in summary_stats.columns and key in ibx_summary.columns]
                        
                        if available_keys:
                            before_merge_cols = len(summary_stats.columns)
                            summary_stats = summary_stats.merge(
                                ibx_summary, 
                                on=available_keys, 
                                how='left'
                            )
                            after_merge_cols = len(summary_stats.columns)
                            print(f"Merge successful: {before_merge_cols} -> {after_merge_cols} columns")
                        else:
                            print(f"WARNING: No common merge keys found for {ibx_category}")
                    else:
                        print(f"WARNING: Empty IBx summary for {ibx_category}")
                        
                except Exception as e:
                    print(f"ERROR processing IBx category {ibx_category}: {e}")
                    import traceback
                    traceback.print_exc()
                    continue              
              
            # Add final summary to collection
            if ibx_category not in all_ibx_dist_dict:
                all_ibx_dist_dict[ibx_category] = {}
            all_ibx_dist_dict[ibx_category][sampling_column] = ibx_dist_dict

        all_summary_dataframes.append(summary_stats)
        print(f"Final summary for {sampling_column}: {summary_stats.shape}")

    if all_inf_ibx:
        all_inf_ibx_df = pd.concat(all_inf_ibx, ignore_index=True)
        all_inf_rh_df  = pd.concat(all_inf_rh, ignore_index=True)
        all_inf_df = pd.merge(all_inf_ibx_df, all_inf_rh_df, on='infIndex', how='outer')
    else:
        all_inf_df = pd.DataFrame()    

    if all_summary_dataframes:
        final_summary = pd.concat(all_summary_dataframes, ignore_index=True)
        print(f"FINAL concatenated summary: {final_summary.shape}")
        print(f"FINAL columns: {list(final_summary.columns)}")
        return final_summary, all_inf_df, all_ibx_dist_dict
    else:
        return pd.DataFrame(), pd.DataFrame(), {}


#####################################################################################
# Heterozygosity calculations
#####################################################################################
def calculate_heterozygosity(maf):
        """Calculate heterozygosity for each locus in the genotype matrix."""
        if maf.ndim == 2:
            maf = np.sum(maf == 1, axis=0) / maf.shape[0]
        else:
            maf = np.array(maf)    

        heterozygosity_func = lambda t: round(1 - (t**2 + (1-t)**2), 4)
        vfunc = np.vectorize(heterozygosity_func)
        site_heterozygosity = vfunc(maf)
        
        return site_heterozygosity


def generate_het_barcode(matrix, indices):
    """Checks for unique alleles at each locus for a specified set of genotypes identified by indices.
    If all alleles are the same at a locus, returns '0' or '1' for that locus.
    If there is a mix of alleles at a locus, returns 'N' for that locus.
    If no indices are provided, returns an empty list.
    
    Note: To update for multi-allelic loci, modify the conditions within the list comprehension.

    TODO: Add option to account for densities to potentially mask polygenomic samples due to low density.
    """
    if isinstance(indices, str):
        indices = ast.literal_eval(indices)

    if len(indices) == 0:
        return 0, [], []
    
    try:
        subset_matrix = matrix[indices, :]
        unique_rows = np.unique(subset_matrix, axis=0)
        
        # Check each column (locus)
        barcode = []
        for col in subset_matrix.T:  # Transpose to iterate over columns
            if np.all(col == 0):
                barcode.append('0')
            elif np.all(col == 1):
                barcode.append('1')
            else:
                barcode.append('N')
        
        het = calculate_heterozygosity(subset_matrix)
        
        return unique_rows.shape[0], barcode, het.tolist()
        
    except Exception as e:
        print(f"Error in generate_het_barcode: {e}")
        return 0, [], []


def process_nested_fws(nested_indices, sampling_df, ibs_matrix = 'ibs_matrix'):
    """Calculate Fws for a single sample based on population heterozygosity.
[        
    Mirrors logic from R package moimix used to calculate F_ws. Specifically following the logic in the function getFws() with this as the comment:

    Compute the within host diversity statistic according to the method devised in  Manske et.al, 2012. Briefly, within sample heterozygosity and within population heterozygosity are computed and assigned to ten equal sized MAF bins [0.0.05]...[0.45,0.5]. For each bin the mean within sample and population heterozygosity is computed. A regression line of these values through the origin is computed for each sample. The \eqn{Fws} is then \eqn{1 - \beta}.
    
    Manske, Magnus, et al. "Analysis of Plasmodium falciparum diversity in natural infections by deep sequencing." Nature 487.7407 (2012): 375-379.
    """
    fws_stats_list = []
    
    # Get year-level indices (same as process_nested_ibx)
    if 'group_year' in nested_indices.keys():
        all_year_indices = nested_indices['group_year']
    elif 'season_bins' in nested_indices.keys():
        all_year_indices = nested_indices['season_bins']
    else:
        print("Warning: No year-level indices found")
        return pd.DataFrame()
    
    # Process each year
    for year_key, year_indices in all_year_indices.items():
        year = str(year_key) if not isinstance(year_key, tuple) else str(year_key[0])
        year_subset = sampling_df[sampling_df['infIndex'].isin(year_indices)]
        
        # Get genome indices for this year
        genome_indices = []
        for idx_list in year_subset['original_nid']:
            if isinstance(idx_list, list):
                genome_indices.extend(idx_list)
        
        if len(genome_indices) == 0:
            continue
            
        # Calculate population-level heterozygosity once per year
        matrix = get_matrix(ibs_matrix)[genome_indices, :]
        group_af = np.sum(matrix == 1, axis=0) / matrix.shape[0]
        group_het = calculate_heterozygosity(matrix)
        maf_bins = pd.cut(group_af, bins=np.linspace(0, 0.5, 11), labels=False) + 1
        group_het_by_bin = pd.Series(group_het).groupby(maf_bins).mean()
        
        # Helper function for Fws calculation
        def calc_fws_for_sample(sample_het_list):
            sample_het = np.array(sample_het_list)
            try:
                sample_het_by_bin = pd.Series(sample_het).groupby(maf_bins).mean()
                combined = pd.DataFrame({
                    'pop_het': group_het_by_bin,
                    'sample_het': sample_het_by_bin
                }).dropna()
                
                if len(combined) == 0:
                    return np.nan
                
                X = combined['pop_het'].values.reshape(-1, 1)
                y = combined['sample_het'].values
                model = LinearRegression(fit_intercept=False)
                model.fit(X, y)
                fws = round(1 - model.coef_[0], 3)
                return fws
            except:
                return np.nan
        
        # Process year-level and nested groups
        for comparison_type, group_data in nested_indices.items():
            if comparison_type in ['group_year', 'season_bins']:
                # Year-level calculation
                year_subset['fws'] = year_subset['heterozygosity'].apply(calc_fws_for_sample)
                valid_fws = year_subset['fws'].dropna()
                
                if len(valid_fws) > 0:
                    fws_summary = _comprehensive_stats(valid_fws, 'fws')
                    fws_stats_list.append({
                        'comparison_type': comparison_type,
                        'year_group': year,
                        'subgroup': None,
                        'allele_frequencies': np.round(group_af, 3).tolist(),
                        'heterozygosity_per_position': np.round(group_het, 3).tolist(),
                        'n_samples': len(valid_fws),
                        **fws_summary.to_dict()
                    })
            else:
                # Handle nested groups within this year
                for key, nested_data in group_data.items():
                    if isinstance(key, tuple):
                        key_year, subgroup = key
                        if str(key_year) == year:
                            if isinstance(nested_data, list) and len(nested_data) > 0:
                                subset_df = year_subset[year_subset['infIndex'].isin(nested_data)]
                                if not subset_df.empty:
                                    # Recalculate group-level stats for this subgroup
                                    subgroup_genome_indices = []
                                    for idx_list in subset_df['original_nid']:
                                        if isinstance(idx_list, list):
                                            subgroup_genome_indices.extend(idx_list)
                                    
                                    if len(subgroup_genome_indices) > 0:
                                        subgroup_matrix = get_matrix(ibs_matrix)[subgroup_genome_indices, :]
                                        subgroup_af = np.sum(subgroup_matrix == 1, axis=0) / subgroup_matrix.shape[0]
                                        subgroup_het = calculate_heterozygosity(subgroup_matrix)
                                        subgroup_maf_bins = pd.cut(subgroup_af, bins=np.linspace(0, 0.5, 11), labels=False) + 1
                                        subgroup_het_by_bin = pd.Series(subgroup_het).groupby(subgroup_maf_bins).mean()
                                        
                                        # Recalculate Fws with subgroup-specific heterozygosity
                                        def calc_fws_subgroup(sample_het_list):
                                            sample_het = np.array(sample_het_list)
                                            try:
                                                sample_het_by_bin = pd.Series(sample_het).groupby(subgroup_maf_bins).mean()
                                                combined = pd.DataFrame({
                                                    'pop_het': subgroup_het_by_bin,
                                                    'sample_het': sample_het_by_bin
                                                }).dropna()
                                                
                                                if len(combined) == 0:
                                                    return np.nan
                                                
                                                X = combined['pop_het'].values.reshape(-1, 1)
                                                y = combined['sample_het'].values
                                                model = LinearRegression(fit_intercept=False)
                                                model.fit(X, y)
                                                return round(1 - model.coef_[0], 3)
                                            except:
                                                return np.nan
                                        
                                        subset_df['fws'] = subset_df['heterozygosity'].apply(calc_fws_subgroup)
                                        valid_fws = subset_df['fws'].dropna()
                                        
                                        if len(valid_fws) > 0:
                                            fws_summary = _comprehensive_stats(valid_fws, 'fws')
                                            fws_stats_list.append({
                                                'comparison_type': comparison_type,
                                                'year_group': year,
                                                'subgroup': str(subgroup),
                                                'allele_frequencies': np.round(subgroup_af, 3).tolist(),
                                                'heterozygosity_per_position': np.round(subgroup_het, 3).tolist(),
                                                'n_samples': len(valid_fws),
                                                **fws_summary.to_dict()
                                            })
    
    if fws_stats_list:
        return pd.DataFrame(fws_stats_list)
    else:
        print("Warning: No Fws summary generated")
        return pd.DataFrame()
    

#####################################################################################
# Matching partner summary statistics
#####################################################################################
# Directly copied from available code for Rh metric paper - not needed for individual sims, but hold here since it's a cross simulation comparison metric. 
# def inverse_var_weight(p_array, var_array):
#     var_weighted_num= np.nansum(np.divide(p_array, 
#                                               var_array, where=var_array!=0))
#     var_weighted_denom = np.nansum(np.divide([1 for _ in var_array], 
#                                                  var_array,where=var_array!=0))
#     weighted_mean = var_weighted_num / var_weighted_denom
    
#     weighted_var = 1/np.sum(1/var_array)
#     weighted_std = np.sqrt(weighted_var)
#     weighted_ci = (weighted_mean - 1.96 * weighted_std,
#                        weighted_mean + 1.96 * weighted_std)
    
    
#     return weighted_mean, weighted_var, weighted_ci


def sample_from_distribution(dist_dict, n_bootstraps = 200, exclude_keys=[1]):
    """ Unpacks the pairwise IBS distribution dictionary to sample from the distribution n times. Excluded keys are values that do not represent the distribution of interest - e.g. IBS=1 for identical barcodes since these would not be detected as a mixed infection. """

    if exclude_keys is None:
        exclude_keys = []
    filtered_dict = {k: v for k, v in dist_dict.items() if k not in exclude_keys}

    values = list(filtered_dict.keys())
    weights = list(filtered_dict.values())

    distribution_list = np.random.choice(values, size=n_bootstraps, p=np.array(weights)/sum(weights))

    return distribution_list


def calculate_individual_rh(barcode_heterozygosity, bootstrap_list):
    """ Calculate individual R_h value based on infection heterozygosity and sampled H_Mono distribution. """
    rh_individual_dist = list(map(lambda i: (i-barcode_heterozygosity)/i if i != 0 else 0, bootstrap_list))
    
    rh_inferred_mean = round(np.median(rh_individual_dist), 3)
    
    return rh_inferred_mean


def calculate_population_rh(df, monogenomic_dict, n_mono_boostraps=200):
    """
    Calculate the R_h statistic for the given sampling dataframe and IBS matrix.

    Args:
        df: DataFrame containing infection data, needs to include effective COI information.
        monogenomic_dict: Dictionary containing the IBS distribution data.
        n_mono_boostraps: Number of bootstrap samples to draw for the population H_Mono estimation.

    Logic from published inferential model (Wong et al 2022, https://doi.org/10.1093/pnasnexus/pgac187):
    - Identify the baseline co-transmission relatedness of unique monogenomic barcode pairs (H_Mono - 200 pairwise draw)
    - Identify per sample polygenomic heterozygosity by the number of Ns (H_Poly)
    - Calculate individual R_h = (H_Mono - H_Poly) / H_Mono
    - Calculate summary statistics (mean, median, std) from the individual polygenomic R_h values per site per year. 

    Adapting for the model:
    - (On hold -tested with too few infections to determine interpretation) For samples super infection samples with an effective COI=2, H_Mono (measured) can be measured directly as defined by the expectation of IBS values between the two unrelated genotypes in a mixed infection.
    - Replicate the bootstrap sampling of H_Mono by drawing from the IBS distribution of monogenomic samples, excluding IBS=1 values (i.e. identical barcodes) to calculate a H_Mono (inferred) value.
    - Polygenomic sample heterozygosity is calculated as the proportion of Ns in the barcode, assuming all alleles in an infection are detectable. Updates to make this more or less sensitive to minor alleles can be made in the generate_het_barcode function.
    """

    
    # coi2_superinfections = df[(df['effective_coi'] == 2) & (df['cotx'] == False)]
    # rh_mono_mean = round(coi2_superinfections['ibs_mean'].mean(), 3)
    poly_samples = df[df['effective_coi'] > 1].copy()
    if 'barcode_N_prop' not in poly_samples.columns:
        poly_samples['barcode_N_prop'] = poly_samples.apply(
            lambda row: row['barcode_with_Ns'].count('N') / len(row['barcode_with_Ns']) if isinstance(row['barcode_with_Ns'], str) else 0, axis=1
        )
    
    bootstrap_list = sample_from_distribution(monogenomic_dict, n_bootstraps=n_mono_boostraps)
    
    poly_samples['individual_inferred_rh'] = poly_samples.apply(lambda row: calculate_individual_rh(row['barcode_N_prop'], bootstrap_list), axis=1)

    rh_measurements = pd.DataFrame([{
        # 'rh_mono_count': len(coi2_superinfections),
        # 'rh_mono_measured_mean': rh_mono_mean,
        # 'rh_mono_measured_median': round(coi2_superinfections['ibs_mean'].median(), 3),
        # 'rh_mono_measured_std': round(coi2_superinfections['ibs_mean'].std(), 3),
        # 'rh_poly_count': len(poly_samples),
        # 'rh_poly_measured_mean': round(poly_samples['individual_measured_rh'].mean(), 3),
        # 'rh_poly_measured_median': round(poly_samples['individual_measured_rh'].median(), 3),
        # 'rh_poly_measured_std': round(poly_samples['individual_measured_rh'].std(), 3),
        'rh_poly_inferred_mean': round(poly_samples['individual_inferred_rh'].mean(), 3),
        'rh_poly_inferred_median': round(poly_samples['individual_inferred_rh'].median(), 3),
        'rh_poly_inferred_std': round(poly_samples['individual_inferred_rh'].std(), 3)
    }])

    return rh_measurements, poly_samples[['infIndex', 'individual_inferred_rh']]


    
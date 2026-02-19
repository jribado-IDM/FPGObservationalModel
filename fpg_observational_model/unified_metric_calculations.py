import pandas as pd
import numpy as np
import ast
from itertools import combinations, chain
from functools import reduce
from collections import Counter
from scipy import stats
from sklearn.linear_model import LinearRegression
import idm

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


def create_empty_structure(source_dict):
    """
    Recursively creates an empty dictionary with the same nested structure 
    as the source_dict.
    """
    empty_dict = {}
    for key, value in source_dict.items():
        if isinstance(value, dict):
            # If the value is a dictionary, recursively call the function
            empty_dict[key] = create_empty_structure(value)
        else:
            empty_dict[key] = None 
    return empty_dict

#####################################################################################
# Group identification for calculations
#####################################################################################
def identify_nested_comparisons(df, time_group, config=None):
    """
    Generate a list of infections within sampling schemes for looping through nested comparisons.
    
    Args:
        df: Input dataframe
        time_group: str — column name defining the time grouping
        config: dict with optional keys: 'polygenomic', 'symptomatic', 'age_bins'
    
    Returns:
        nested_indices[population][sampling_column][grouping] = {time_key: [infIndex, ...]}
    """

    grouping_dict = {}

    # --- All infections by time group ---
    grouping_dict['all'] = df.groupby(time_group)['infIndex'].apply(list).to_dict()

    # --- Polygenomic ---
    if config.get('polygenomic', False):
        subset = df.copy()
        subset['is_polygenomic'] = subset['effective_coi'] > 1
        poly_vals = subset['is_polygenomic'].unique()
        if True in poly_vals and False in poly_vals:
            grouping_dict['polygenomic'] = (
                subset.groupby([time_group, 'is_polygenomic'])['infIndex']
                .apply(list).to_dict()
            )
        else:
            print(f"User specified polygenomic comparisons, but only one group available.")

    # --- Symptomatic ---
    if config.get('symptomatic', False):
        if len(df['fever_status'].unique()) > 1:
            grouping_dict['symptomatic'] = (
                df.groupby([time_group, 'fever_status'])['infIndex']
                .apply(list).to_dict()
            )
        else:
            print(f"User specified fever status comparisons, but only one status available.")

    # --- Age bins from config ---
    if config.get('age_bins', False):
        days_per_year = 365.25
        subset = df.copy()
        age_bins = [0, int(days_per_year * 5), int(days_per_year * 15), int(subset['age_day'].max() + 1)]
        age_bin_labels = ['0-5yrs', '5-15yrs', '15+yrs']
        subset['age_bin'] = pd.cut(
            subset['age_day'], bins=age_bins, labels=age_bin_labels, include_lowest=True
        )
        if len(subset['age_bin'].unique()) > 1:
            grouping_dict['age_bins'] = (
                subset.groupby([time_group, 'age_bin'], observed=True)['infIndex']
                .apply(list).to_dict()
            )
        else:
            available = subset['age_bin'].unique()[0]
            print(f"User specified age bin comparisons, but only one age group available: {available}.")

    return grouping_dict

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
    # Add length check
    if len(series) <= 1:
        val= series.iloc[0] if len(series) == 1 else np.nan
        return pd.Series({
            f'{prefix}_mean': val,
            f'{prefix}_median': val,
            f'{prefix}_std': 0.0,
            f'{prefix}_q25': val,
            f'{prefix}_q75': val,
            f'{prefix}_min': val,
            f'{prefix}_max': val
        })
    
    stats = series.describe()
    result = {
        f'{prefix}_mean': round(stats['mean'], 3),
        f'{prefix}_median': round(stats['50%'], 3),
        f'{prefix}_std': round(stats['std'], 3),
        f'{prefix}_q25': round(stats['25%'], 3),
        f'{prefix}_q75': round(stats['75%'], 3),
        f'{prefix}_min': round(stats['min'], 3),
        f'{prefix}_max': round(stats['max'], 3)
    }

    if any(x in prefix for x in ['coi', 'true', 'effective', 'genotype']):
        poly_mask = series > 1
        count = poly_mask.sum()
        result[f'{prefix}_poly_count'] = count
        result[f'{prefix}_poly_prop'] = round(count / len(series), 3)

    return pd.Series(result)


def _analyze_binary_in_subset(series, prefix):
    """Analyze binary column within a subset defined by mask."""
    if len(series) == 0:
        return {f'{prefix}_count': 0, f'{prefix}_prop': np.nan}
    
    mask = series == 1 
    count = mask.sum()
    prop = round(count / len(series), 3)
    
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
    
    if count > 1:
        summary_data = {
            f'{ibx_prefix}_count': int(count),
            f'{ibx_prefix}_mean': round(np.mean(expanded_values), 3),
            f'{ibx_prefix}_std': round(np.std(expanded_values, ddof=1), 3), # Use ddof=1 for sample std to match pandas
            f'{ibx_prefix}_q25': round(np.percentile(expanded_values, 25), 3),
            f'{ibx_prefix}_median': round(np.median(expanded_values), 3),
            f'{ibx_prefix}_q75': round(np.percentile(expanded_values, 75), 3),
            f'{ibx_prefix}_min': round(np.min(expanded_values), 3),
            f'{ibx_prefix}_max': round(np.max(expanded_values), 3)
        }
    else:
        val = expanded_values[0]
        summary_data = {        
            f'{ibx_prefix}_count': int(count),
            f'{ibx_prefix}_mean': val,
            f'{ibx_prefix}_std': 0.0, 
            f'{ibx_prefix}_q25': val,
            f'{ibx_prefix}_median': val,
            f'{ibx_prefix}_q75': val,
            f'{ibx_prefix}_min': val,
            f'{ibx_prefix}_max': val
    }
    return pd.DataFrame([summary_data]) 


#####################################################################################
# Run summaries and IBx calculations for nested groups
#####################################################################################
def process_nested_summaries(nested_indices, sampling_df):
   
    summary_stats_list = []

    for grouping, data in nested_indices.items():
        for key, indices in data.items():
            if isinstance(key, tuple):
                time_key = key[0]
                subpopulation_group = key[1]
            else:
                time_key = key
                subpopulation_group = None    

            group_subset = sampling_df[sampling_df['infIndex'].isin(indices)]
            summary = comprehensive_group_summary(group_subset)
            summary = summary.to_dict() if isinstance(summary, pd.Series) else summary

            summary_dict = {
                'time_value': time_key,
                'comparison_type': grouping,
                'comparison_group': subpopulation_group,
                **summary
            }            
            summary_stats_list.append(summary_dict)
    
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
    save_ibx_distributions=True,
    save_pairwise_ibx=False):
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

    ibx_dist_dict, individual_ibx_dict = {}, {}
    ibx_summ_list = []
    ibx_results_df, individual_ibx_df = pd.DataFrame(), pd.DataFrame()

    all_year_indices = nested_indices['all']    
    for time_key, indices in all_year_indices.items():     
        time_subset = df[df['infIndex'].isin(indices)]
        time_subset = update_ibx_index(time_subset)

        # Step 1: Run pairwise IBx calculations once per year
        genome_indices = []
        for idx_list in time_subset['recursive_nid']:
            if isinstance(idx_list, list):
                genome_indices.extend(idx_list)

        matrix = get_matrix(gt_matrix)[genome_indices, :]
        print("Genotype matrix shape:", matrix.shape)
        ibx_indices, ibx_matrix = calculate_ibx_matrix(time_subset, matrix)

        # FOR HAIRBALL connectedness plots, save the ibx_matrix and ibx_indices per year 
        if save_pairwise_ibx:
            # Update output directory here
            output_dir = "output"
            pd.save_csv(ibx_indices, f"{output_dir}/ibx_indices_{time_key}.csv", index=False)
            np.save(f"{output_dir}/ibx_matrix_{time_key}.npy", ibx_matrix)

        # Add column with the ibx_index for each infection
        ibx_mapping = dict(zip(ibx_indices['ibx_nid'], ibx_indices['ibx_index']))
        time_subset['ibx_index'] = time_subset['ibx_nid'].apply(lambda nid_list: [ibx_mapping[nid] for nid in nid_list] if isinstance(nid_list, list) else None)
        
        # Step 2: Run IBx summaries for each year
        for grouping, data in nested_indices.items():
            if grouping not in ibx_dist_dict:  # only initialize if it doesn't exist yet
                ibx_dist_dict[grouping] = {}
            for key, values in data.items():
                if isinstance(key, tuple):
                    if time_key == key[0]:
                        subpopulation_group = key[1]
                        subset_infections = values
                    else:
                        continue  
                else:
                    if key != time_key:
                        continue  
                    subpopulation_group = None
                    subset_infections = values  

                subset_df = time_subset[time_subset['infIndex'].isin(subset_infections)]
                subset_indices = list(chain.from_iterable(subset_df['ibx_index'].tolist()))         

                if isinstance(subset_indices, list) and len(subset_indices) > 1:
                    distribution = ibx_distribution(subset_indices, ibx_matrix)
                    summary_stats = weighted_describe_scipy(distribution, f"pop-{ibx_prefix}")
                    
                    result_row = {
                        'time_value': time_key,
                        'comparison_type': grouping,
                        'comparison_group': subpopulation_group,
                    }
                    result_row.update(summary_stats.iloc[0].to_dict())
                    ibx_summ_list.append(result_row)

                    if save_ibx_distributions:
                        ibx_dist_dict[grouping][key] = distribution
        
            
        # Step 3: Individual IBx calculations for polygenomic infections
        if individual_ibx_calculation:
            polygenomic_subset = time_subset[time_subset['effective_coi'] > 1]
            if not polygenomic_subset.empty:
                polygenomic_dict = dict(zip(polygenomic_subset['infIndex'], polygenomic_subset['ibx_index']))

                for inf_id, ibx_list in polygenomic_dict.items():
                    distribution = ibx_distribution(ibx_list, ibx_matrix)
                    individual_ibx_dict[inf_id] = weighted_describe_scipy(distribution, ibx_prefix) 

    if ibx_summ_list:
        ibx_results_df = pd.DataFrame(ibx_summ_list)

    if individual_ibx_dict:  
        individual_ibx_df = pd.concat(individual_ibx_dict, names=['infIndex', 'row_id']).reset_index(level=0)

    return ibx_results_df, individual_ibx_df, ibx_dist_dict


def process_individual_ibx(nested_dict, individual_ibx_df, ibx_prefix):
    """ Merge individual IBx calculations back to the sampling dataframe based on nested comparison groups. """
    if individual_ibx_df.empty:
        print("No individual IBx data to process.")
        return pd.DataFrame()
    
    merged_df = pd.DataFrame()
    for grouping, data in nested_dict.items():
        if grouping == "polygenomic":
            continue
        for key, indices in data.items():
            if isinstance(key, tuple):
                time_key = key[0]
                subpopulation_group = key[1]
            else:
                time_key = key
                subpopulation_group = None    

            time_subset = individual_ibx_df[individual_ibx_df['infIndex'].isin(indices)]
            
            result = _comprehensive_stats(time_subset[f"{ibx_prefix}_mean"], f"ind-{ibx_prefix}").to_dict()
            summary = {
                'time_value': time_key,
                'comparison_type': grouping,
                'comparison_group': subpopulation_group,
                **result
                }
    
            merged_df = pd.concat([merged_df, pd.DataFrame([summary])], ignore_index=True)            
           
    return merged_df


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

    # Helper function to generate the heterozygosity barcode for a sample
    # Calculate population-level heterozygosity once per year
    def calc_group_maf(genome_indices, ibs_matrix = 'ibs_matrix'): 
        matrix = get_matrix(ibs_matrix)[genome_indices, :]
        group_af = np.sum(matrix == 1, axis=0) / matrix.shape[0]
        group_maf = np.minimum(group_af, 1 - group_af) 
        group_het = calculate_heterozygosity(matrix)
        maf_bins = pd.cut(group_maf, bins=np.linspace(0, 0.5, 11), labels=False) + 1
        group_het_by_bin = pd.Series(group_het).groupby(maf_bins).mean()

        return group_af, group_het, maf_bins, group_het_by_bin

    # Helper function for Fws calculation
    def calc_fws_for_sample(sample_het_list, maf_bins, het_bins):
        sample_het = np.array(sample_het_list)
        try:
            sample_het_by_bin = pd.Series(sample_het).groupby(maf_bins).mean()
            combined = pd.DataFrame({
                'pop_het': het_bins,
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


    for grouping, data in nested_indices.items():
        for key, indices in data.items():
            if isinstance(key, tuple):
                time_key = key[0]
                subpopulation_group = key[1]
            else:
                time_key = key
                subpopulation_group = None    

            time_subset = sampling_df[sampling_df['infIndex'].isin(indices)]

            genome_indices = []
            for idx_list in time_subset['original_nid']:
                if isinstance(idx_list, list):
                    genome_indices.extend(idx_list)

            if len(genome_indices) == 0:
                continue

            group_af, group_het, maf_bins, group_het_by_bin = calc_group_maf(genome_indices, ibs_matrix)

            # Time-level calculation
            time_subset = time_subset.copy()
            time_subset['fws'] = time_subset['heterozygosity'].apply(lambda x: calc_fws_for_sample(x, maf_bins, group_het_by_bin))
                        
            fws_stats_dict = {
                'time_value': time_key,
                'comparison_type': grouping,
                'comparison_group': subpopulation_group,
                'allele_frequencies': np.round(group_af, 2).tolist(),
                'heterozygosity_per_position': np.round(group_het, 2).tolist()
                } 

            valid_fws = time_subset['fws'].dropna()
            if len(valid_fws) > 0:
                fws_summary = _comprehensive_stats(valid_fws, 'fws')
                fws_stats_list.append({
                    **fws_stats_dict,
                    **fws_summary.to_dict()
                })
            else:
                fws_stats_list.append({**fws_stats_dict,
                    'fws_mean': 0.0,
                    'fws_median': 0.0,
                    'fws_std': 0.0,
                    'fws_q25': 0.0,
                    'fws_q75': 0.0,
                    'fws_min': 0.0,
                    'fws_max': 0.0
                })    

    if fws_stats_list:
        return pd.DataFrame(fws_stats_list)
    else:
        print("Warning: No Fws summary generated")
        return pd.DataFrame()
    

#####################################################################################
# Matching partner summary statistics
#####################################################################################
def sample_from_distribution(dist_dict, n_bootstraps = 200, exclude_keys=[1]):
    """ Unpacks the pairwise IBS distribution dictionary to sample from the distribution n times. Excluded keys are values that do not represent the distribution of interest - e.g. IBS=1 for identical barcodes since these would not be detected as a mixed infection. """

    if exclude_keys is None:
        exclude_keys = []
    filtered_dict = {k: v for k, v in dist_dict.items() if k not in exclude_keys}

    values = list(filtered_dict.keys())
    weights = list(filtered_dict.values())

    if filtered_dict:
        distribution_list = np.random.choice(values, size=n_bootstraps, p=np.array(weights)/sum(weights))
    else:
        # Account for extreme scenarios where there are no values to sample from after filtering (e.g. all IBS=1 pairs, clonal population)
        distribution_list = np.array([])    

    return distribution_list


def calculate_individual_rh(barcode_heterozygosity, bootstrap_list):
    """ Calculate individual R_h value based on infection heterozygosity and sampled H_Mono distribution. """
    rh_individual_dist = list(map(lambda i: (i-barcode_heterozygosity)/i if i != 0 else 0, bootstrap_list))
    
    rh_inferred_mean = round(np.mean(rh_individual_dist), 3)
    
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

    poly_samples = df[df['effective_coi'] > 1].copy()
    if 'barcode_N_prop' not in poly_samples.columns:
        poly_samples['barcode_N_prop'] = poly_samples.apply(
            lambda row: row['barcode_with_Ns'].count('N') / len(row['barcode_with_Ns']) if isinstance(row['barcode_with_Ns'], str) else 0, axis=1
        )
    
    bootstrap_list = sample_from_distribution(monogenomic_dict, n_bootstraps=n_mono_boostraps)
    
    if len(bootstrap_list) > 0:
        poly_samples['individual_inferred_rh'] = poly_samples.apply(lambda row: calculate_individual_rh(row['barcode_N_prop'], bootstrap_list), axis=1)
        rh_measurements = {
        'rh_poly_inferred_mean': round(poly_samples['individual_inferred_rh'].mean(), 3),
        'rh_poly_inferred_median': round(poly_samples['individual_inferred_rh'].median(), 3),
        'rh_poly_inferred_std': round(poly_samples['individual_inferred_rh'].std(), 3)
        }

    else:    
        poly_samples['individual_inferred_rh'] = np.nan
        rh_measurements = {
            'rh_poly_inferred_mean': np.nan,
            'rh_poly_inferred_median': np.nan,
            'rh_poly_inferred_std': np.nan
        }

    return rh_measurements, poly_samples[['infIndex', 'individual_inferred_rh']]


def process_nested_rh(nested_dict, sampling_df, ibx_dist_dict, inf_ibx):
    """
    Calculate R_h for nested comparison groups.
    """

    yearly_rh_df, inf_rh = [], pd.DataFrame()

    # Get all the distribution dictionaries for H_mono bootstrap calculations
    monogenomic_dict = {k[0]: v for k,v in ibx_dist_dict.items() if k[1]==False}
    
    for time_key, indices in nested_dict['all'].items():
        if time_key not in monogenomic_dict.keys():
            print(f"Warning: No monogenomic IBS distribution found for time group {time_key}. R_h calculations will be skipped for this group.")
            continue
        else:
            df = sampling_df[sampling_df['infIndex'].isin(indices)]  
            df = df.merge(pd.DataFrame(inf_ibx), on='infIndex', how='left')

            rh_summary, sample_rh = calculate_population_rh(df, monogenomic_dict[time_key])

            rh_combined_summary = {
                'time_value': time_key,
                'comparison_type': 'all',
                'comparison_group': None,
                **rh_summary
            }
            yearly_rh_df.append(rh_combined_summary)
            inf_rh = pd.concat([inf_rh, sample_rh], ignore_index=True)

    return pd.DataFrame(yearly_rh_df), inf_rh
    


################################################################################# Putting it all together
################################################################################
def run_time_summaries(sample_df,
                        subpop_config=None,
                        user_ibx_categories=None,
                        individual_ibx_calculation=True,
                        fws_calculation=True,
                        rh_calculation=True,
                        save_ibx_distributions=True):

    df = sample_df
    sampling_columns = df.filter(regex="rep").columns.to_list()

    all_summary_dataframes, all_inf_ibx_chunks = [], []
    all_ibx_dist_dict = {ibx_category: {} for ibx_category in user_ibx_categories} if user_ibx_categories else {}

    for sampling_column_name in sampling_columns:
        sampling_df = df[df[sampling_column_name].notna()]
        print(sampling_df.head())
        print(f"Starting: {sampling_column_name}")

        # --- Resolve time group ---
        if 'month' in sampling_column_name:
            time_group = 'group_month'
        elif 'random' in sampling_column_name:
            time_group = 'group_year'
        elif 'seasonal' in sampling_column_name:
            if len(df[sampling_column_name].unique()) > 1:
                sampling_df['group_season'] = sampling_df[sampling_column_name]
                time_group = 'group_season'
            else:
                print(f"User specified comparisons by season, but only one season found.")
        else:
            time_group = sampling_column_name
            print(f"User specified comparisons by other grouping column, {sampling_column_name}.")

        for pop_key, population_subset in sampling_df.groupby('population'):

            nested_dict = identify_nested_comparisons(population_subset, time_group, config=subpop_config)

            summary_stats = process_nested_summaries(nested_dict, population_subset)
            summary_stats.insert(0, 'population', pop_key)
            summary_stats.insert(1, 'sampling_scheme', sampling_column_name)
            summary_stats.insert(2, 'time_group', time_group)

            if fws_calculation:
                fws_stats_timescale = process_nested_fws(nested_dict, sampling_df)
                if not fws_stats_timescale.empty:
                    summary_stats = summary_stats.merge(fws_stats_timescale,
                                                        on=['time_value', 'comparison_type', 'comparison_group'],
                                                        how='left')
                else:
                    print(f"Warning: No Fws stats generated for {sampling_column_name} for population {pop_key}")

            if time_group == 'group_month':
                print(f"IBx calculations are only run for year-level comparisons. Skipping for {time_group}.")
            else:
                if user_ibx_categories and len(user_ibx_categories) > 0:

                    sample_inf_rh = []
                    sample_inf_ibx = {}

                    for ibx_category in user_ibx_categories:
                        ibx_dist_dict = {}
                        print(f"\nProcessing IBx category: {ibx_category} | population: {pop_key} | time group: {time_group}")

                        try:
                            ibx_summary, ibx_inf, ibx_dist_dict = process_nested_ibx(
                                population_subset,
                                f'{ibx_category}_matrix',
                                nested_dict,
                                ibx_prefix=ibx_category,
                                individual_ibx_calculation=individual_ibx_calculation,
                                save_ibx_distributions=save_ibx_distributions
                            )

                            if not ibx_inf.empty:
                                within_inf_summary = process_individual_ibx(nested_dict, ibx_inf, ibx_category)
                                ibx_summary = ibx_summary.merge(within_inf_summary,
                                    on=['time_value', 'comparison_type', 'comparison_group'],
                                    how='left')
                                if ibx_category in sample_inf_ibx and not sample_inf_ibx[ibx_category].empty:
                                    sample_inf_ibx[ibx_category] = pd.concat([sample_inf_ibx[ibx_category], ibx_inf], ignore_index=True)
                                else:
                                    sample_inf_ibx[ibx_category] = ibx_inf

                            summary_stats = summary_stats.merge(ibx_summary,
                                on=['time_value', 'comparison_type', 'comparison_group'],
                                how='left')

                            if rh_calculation and ibx_category == 'ibs':
                                if 'polygenomic' not in nested_dict.keys():
                                    print(f"Warning: No polygenomic subpopulation comparisons found for {pop_key}.")
                                    if len(population_subset['effective_coi'].unique()) != 1:
                                        print(f"Warning: No polygenomic infections found for {pop_key}.")
                                        continue
                                else:
                                    rh_summary, sample_rh = process_nested_rh(nested_dict, population_subset, ibx_dist_dict['polygenomic'], ibx_inf)
                                    summary_stats = summary_stats.merge(rh_summary, on=['time_value', 'comparison_type', 'comparison_group'], how='left')
                                    sample_inf_rh.append(sample_rh)

                        except Exception as e:
                            print(f"ERROR processing IBx category {ibx_category} / {pop_key}: {e}")
                            import traceback
                            traceback.print_exc()
                            continue

                        if ibx_dist_dict:
                            if f"population_{pop_key}" not in all_ibx_dist_dict[ibx_category]:
                                all_ibx_dist_dict[ibx_category][f"population_{pop_key}"] = {}
                            all_ibx_dist_dict[ibx_category][f"population_{pop_key}"][sampling_column_name] = ibx_dist_dict
                    
                    if sample_inf_ibx:
                        dfs = list(sample_inf_ibx.values())
                        if dfs:
                            merged = reduce(lambda left, right: left.merge(right, on='infIndex', how='outer'), dfs)
                            del dfs  # free memory
                            
                            if sample_inf_rh:
                                merged = merged.merge(pd.concat(sample_inf_rh, ignore_index=True), on='infIndex', how='left')
                            all_inf_ibx_chunks.append(merged)
                            del merged  # free the local reference

            all_summary_dataframes.append(summary_stats)
            print(f"Final summary for {sampling_column_name} / {pop_key}: {summary_stats.shape}")

    # Concat all IBx chunks once at the end
    if all_inf_ibx_chunks:
        all_inf_ibx_df = pd.concat(all_inf_ibx_chunks, ignore_index=True)
    else:
        all_inf_ibx_df = pd.DataFrame()

    if all_summary_dataframes:
        final_summary = pd.concat(all_summary_dataframes, ignore_index=True)
        print(f"FINAL concatenated summary: {final_summary.shape}")
        print(f"FINAL columns: {list(final_summary.columns)}")
        return final_summary, all_inf_ibx_df, all_ibx_dist_dict
    else:
        return pd.DataFrame(), pd.DataFrame(), {}
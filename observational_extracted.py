import pandas as pd
import numpy as np
import ast
import math
from concurrent.futures import ProcessPoolExecutor, as_completed

def parse_list(s):
    """
    Converts a string representation of a list into an actual Python list.
    If conversion fails, returns an empty list.
    """
    try:
        return ast.literal_eval(s)
    except Exception:
        return []

def assign_season_group(row):
    year = int(row['year'])
    month = int(row['month'])
    if month == 1:
        # For January, group with the previous year's August-December period.
        return f"Wet season: {year-1}-08 to {year}-01"
    elif month >= 8:
        # For months Aug-Dec, group with January of the following year.
        return f"Wet season: {year}-08 to {year+1}-01"
    else:
        # For other months (e.g., February to July), you could leave it ungrouped.
        return f"Dry season: {year}-02 to {year}-07"

def assign_peak_group(row):
    year = int(row['year'])
    month = int(row['month'])
    if month >= 10 and month <= 12:
        return f"Wet season: {year}-10 to {year}-12"
    elif 3 <= month <= 6:
        return f"Dry season: {year}-03 to {year}-06"
    else:
        return f"Off-peak season"      


def reassign_by_intervention(df, intervention_reset=29):
    df['intervention_month'] = df['continuous_month'] - intervention_reset
    df['year'] = df['intervention_month'] // 12
    return(df)


input_seed=418
def process_genetic_data(df):
    """
    Processes the DataFrame by:
      1. Parsing the 'genome_ids' column and computing:
         - true_coi: the total count of items in genome_ids.
         - effective_coi: the count of unique items in genome_ids.
      2. Parsing the 'bite_ids' column and computing:
         - cotransmission: a Boolean indicating if all items in bite_ids are unique.
    
    Parameters:
      df (pd.DataFrame): Input DataFrame containing at least 'genome_ids' and 'bite_ids' columns as strings.
    
    Returns:
      pd.DataFrame: The modified DataFrame with additional computed columns.
    """

    # Sample one representative infection per person per year
    df = df.groupby(['IndividualID', 'simulation_year']).sample(n=1, random_state=input_seed)
    
    # Sample one representative infection per person per month
    #df = df.groupby(['IndividualID', 'month']).sample(n=1, random_state=input_seed)

    # Sample peak and trough timepoints
    # peaks = [668, 1033, 1398, 1763, 2128] 
    # troughs = [501, 866, 1231, 1596, 1961]
    # df = df[df['day'].isin(peaks + troughs)]

    # 1. Process the 'genome_ids' column
    df["genome_ids"] = df["genome_ids"].apply(parse_list)
    df["true_coi"] = df["genome_ids"].apply(len)
    df["effective_coi"] = df["genome_ids"].apply(lambda x: len(set(x)))
    
    # 2. Process the 'bite_ids' column 
    df["bite_ids"] = df["bite_ids"].apply(parse_list)
    df["superinfection"] = df["bite_ids"].apply(lambda x: len(set(x)) > 1)
    df["cotransmission"] = df["bite_ids"].apply(lambda x: len(set(x)) == 1)

    df["season"] = df.apply(assign_season_group, axis=1)
    df["peak_season"] = df.apply(assign_peak_group, axis=1)
    
    return df


def summarize_infections(df, groupby_cols=['year', 'month'], 
    sample_n=None, sample_proportionally=True, sample_seasons=False,
    seed=input_seed):
    """
    Group and summarize infection data by a given time frame.
    
    For each group (e.g., by year and month), computes:
      1. Total rows, count and proportion of rows with effective_coi > 1.
      2. Total and unique counts (and the proportion) of genome_ids 
         (after flattening the lists from all rows).
      3. Count and proportion for rows where cotransmission and superinfection is True.
    
    Optionally, a random sample of rows can be taken before grouping.
    
    Parameters:
      df (pd.DataFrame): Input DataFrame. Must contain columns:
                         - 'effective_coi' (numeric, computed beforehand)
                         - 'genome_ids' (list-like; if not, parse first with parse_list)
                         - 'cotransmission' (Boolean)
                         - Time frame columns (default: 'year' and 'month')
      groupby_cols (list): List of column names to group by. Defaults to ['year', 'month'].
      sample_n (int or None): If provided, randomly sample n rows from the DataFrame before summarizing.
      seed (int, optional): Random seed for sampling.
    
    """
    # Optionally, subsample the DataFrame.
    if sample_n is not None:
        if sample_proportionally:
            if not sample_seasons:
                # Group by 'year'; each group will sample up to sample_n rows, but if the group
                # has fewer than sample_n rows, it will take all available rows.
                df = df.groupby(['year'], group_keys=False).apply(lambda grp: grp.sample(n=min(len(grp), sample_n), random_state=seed))
            else:
                sample_season = math.floor(sample_n)
                df = df.groupby(['season'], group_keys=False).apply(lambda grp: grp.sample(n=min(len(grp), sample_season), random_state=input_seed))
        else:
            # Group by both 'year' and 'month'
            sample_monthly = math.floor(sample_n/12)
            df = df.groupby(['continuous_month'], group_keys=False).apply(lambda grp: grp.sample(n=min(len(grp), sample_monthly), random_state=input_seed))
    
    def group_summary(group):
        n = len(group)
        # Count poly_coi: rows with effective_coi > 1.
        poly_count = (group['effective_coi'] > 1).sum()
        poly_prop = poly_count / n if n > 0 else None

        # effective coi mean
        effective_coi_mean = group['effective_coi'].mean()
        true_coi_mean = group['true_coi'].mean()
        
        # Flatten the genome_ids lists from all rows in the group.
        all_genome_ids = [gid for sublist in group['genome_ids'] 
                                 if isinstance(sublist, list)
                                 for gid in sublist]
        total_genome = len(all_genome_ids)
        unique_genome = len(set(all_genome_ids))
        unique_prop = unique_genome / total_genome if total_genome > 0 else None
        
        # Count cotransmission True rows.
        polygenomic = group[group['effective_coi'] > 1]
        cotrans_count = polygenomic['cotransmission'].sum()  # Assuming boolean where True==1, False==0
        cotrans_prop = cotrans_count / poly_count if poly_count > 0 else None

        supertrans_count = polygenomic['superinfection'].sum()
        supertrans_prop = supertrans_count / poly_count if poly_count > 0 else None
        
        return pd.Series({
            'n_infections': n,
            'poly_coi_count': poly_count,
            'poly_coi_prop': poly_prop,
            'true_coi_mean': true_coi_mean,
            'effective_coi_mean': effective_coi_mean,
            'genome_ids_total_count': total_genome,
            'genome_ids_unique_count': unique_genome,
            'genome_ids_unique_prop': unique_prop,
            'cotransmission_count': cotrans_count,
            'cotransmission_prop': cotrans_prop,
            'superinfection_count': supertrans_count,
            'superinfection_prop': supertrans_prop
        })
    
    summary_df = df.groupby(groupby_cols).apply(group_summary).reset_index()
    return summary_df

def combined_summaries(df):
    """
    Generate combined summaries of infection data using different sampling schemes.
    
    Parameters:
      df (pd.DataFrame): Input DataFrame. Must contain columns:
                         - 'effective_coi' (numeric, computed beforehand)
                         - 'genome_ids' (list-like; if not, parse first with parse_list)
                         - 'cotransmission' (Boolean)
                         - Time frame columns (default: 'year' and 'month')
    
    Returns:
      pd.DataFrame: A combined summary DataFrame with different sampling schemes.
    """
    # Summarize using different sampling schemes
    input_n_sample = 100
    summary_all_yearly=summarize_infections(df, groupby_cols=['year']).assign(sampling_scheme='All - Yearly', sample_n=None)
    summary_all_monthly=summarize_infections(df, groupby_cols=['continuous_month'], sample_n=None).assign(sampling_scheme='All - Monthly')
    #summary_yearly_proportionally = summarize_infections(df, groupby_cols=['year'], sample_n=input_n_sample).assign(sampling_scheme='Sample - Proportional')
    #summary_yearly_evenly = summarize_infections(df, groupby_cols=['year'], sample_n=input_n_sample, sample_proportionally=False).assign(sampling_scheme='Sample - Even')
    summary_yearly_seasonally = summarize_infections(df, groupby_cols=['season'], sample_n=input_n_sample, sample_seasons=True).assign(sampling_scheme='Sample - Seasonal')
    #summary_yearly_peaks = summarize_infections(df, groupby_cols=['peak_season'], sample_n=input_n_sample).assign(sampling_scheme='Sample - Peak Seasonal')

    # combined_df = pd.concat([summary_all_yearly, summary_all_monthly, \
    #    summary_yearly_proportionally, summary_yearly_evenly, \
    #    summary_yearly_seasonally, summary_yearly_peaks], 
    # axis=0)

    combined_df = pd.concat([summary_all_yearly, summary_all_monthly, summary_yearly_seasonally], axis=0)
    #combined_df = summary_all_monthly
    
    return(combined_df)


def process_file(row, output_summary_dir,
reassign_intervention_time=True):
    """Process a single file and write the summary output.
    
    Parameters:
      row (pd.Series): A row from the file list DataFrame.
      output_summary_dir (str): Folder where output files will be saved.
    
    Returns:
      str: The path to the written summary file.
    """
    output_name = row['output_name']
    # Construct the full path to the input file.
    input_file = os.path.join(row['input_dir'], "infIndexRecursive-genomes-df.csv")
    
    # Read and process the input file.
    df = pd.read_csv(input_file)
    df['simulation_year'] = df['year'].copy()
    df['month'] = df['month'] + 1
    df['continuous_month'] = (df["simulation_year"]) * 12 + df["month"]

    if reassign_intervention_time:
        df = reassign_by_intervention(df)
     
    df = process_genetic_data(df)
    run_summary = combined_summaries(df)
    
    # Construct output file path.
    output_file = os.path.join(output_summary_dir, f"{output_name}_summary.csv")
    run_summary.to_csv(output_file, index=False)
    
    return output_file


def process_all_files(file_list_path, output_summary_dir, max_workers=None):
    # Read the file list CSV.
    file_list_df = pd.read_csv(file_list_path)
    
    # Ensure the output directory exists.
    os.makedirs(output_summary_dir, exist_ok=True)

    # Create a list of rows to process.
    # Using DataFrame.iterrows() returns each row as a Series.
    file_rows = [row for _, row in file_list_df.iterrows()]

    # Use ProcessPoolExecutor to parallelize file processing.
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        future_to_file = {executor.submit(process_file, row, output_summary_dir): row['output_name']
                          for row in file_rows}#[:1]}
        
        # Optionally, process the results as they complete.
        for future in as_completed(future_to_file):
            try:
                result = future.result()
                print(f"Processed: {result}")
            except Exception as exc:
                output_name = future_to_file[future]
                print(f"{output_name} generated an exception: {exc}")


project_dir = "/mnt/data/malaria/synthetic_genomes/jessica_projects/2504_GRSweep"
file_list_path = os.path.join(project_dir, "sim_mapping_new_itns_6yr_highBiting.csv")   
# A folder for saving output summary files.
output_summary_dir = os.path.join(project_dir, "infectionFPGReport_summaries_6yr_highBiting")  # Change to your desired output directory.
os.makedirs(output_summary_dir, exist_ok=True)
process_all_files(file_list_path, output_summary_dir, max_workers=4)



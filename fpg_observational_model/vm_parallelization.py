#####################################################################################
# Functions to parallelize the processing of multiple FPG simulations through the observational model on a VM. The goal is to replace this with a dtk_post_process script in the future. 
#####################################################################################
import pandas as pd
import os
from run_observational_model import process_file
from concurrent.futures import ProcessPoolExecutor, as_completed
import logging
from datetime import datetime

# Set up logging
def setup_logging(output_dir):
    """Set up logging for parallel processing."""
    # Ensure output directory exists before creating log file
    os.makedirs(output_dir, exist_ok=True)
    
    log_file = os.path.join(output_dir, f"parallel_processing_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()
        ]
    )
    return logging.getLogger(__name__)

def check_output_exists(sim_name, output_summary_dir, required_files=None):
    """
    Check if output files already exist for a simulation.
    
    Parameters:
        sim_name (str): Name of the simulation
        output_summary_dir (str): Base output directory
        required_files (list): List of required file suffixes to check for
        
    Returns:
        bool: True if all required outputs exist, False otherwise
    """
    if required_files is None:
        required_files = ['_FPG_ModelSummaries.csv', '_FPG_SampledInfections.csv']
    
    sim_output_dir = os.path.join(output_summary_dir, sim_name)
    
    if not os.path.exists(sim_output_dir):
        return False
    
    for suffix in required_files:
        expected_file = os.path.join(sim_output_dir, f"{sim_name}{suffix}")
        if not os.path.exists(expected_file):
            return False
    
    return True

def filter_processed_files(file_list_df, output_summary_dir, check_outputs=True):
    """
    Filter out files that have already been processed.
    
    Parameters:
        file_list_df (pd.DataFrame): DataFrame with simulation info
        output_summary_dir (str): Base output directory
        check_outputs (bool): Whether to check for existing output files
        
    Returns:
        tuple: (to_process_df, already_processed_df)
    """
    if not check_outputs:
        return file_list_df, pd.DataFrame()
    
    to_process = []
    already_processed = []
    
    for _, row in file_list_df.iterrows():
        sim_name = row['output_name']
        
        if check_output_exists(sim_name, output_summary_dir):
            already_processed.append(row)
        else:
            to_process.append(row)
    
    to_process_df = pd.DataFrame(to_process) if to_process else pd.DataFrame()
    already_processed_df = pd.DataFrame(already_processed) if already_processed else pd.DataFrame()
    
    return to_process_df, already_processed_df

def process_all_files(file_list_path, output_summary_dir, config_path=None, max_workers=None, verbose=False, skip_existing=True):
    """
    Process all files in parallel using multiprocessing.
    
    Parameters:
        file_list_path (str): Path to CSV file with 'output_name' and 'input_dir' columns
        output_summary_dir (str): Directory to save all outputs
        config_path (str): Path to config file (optional, will use defaults if not provided)
        max_workers (int): Maximum number of parallel processes (None = use all CPUs)
        verbose (bool): Whether to enable verbose logging
        skip_existing (bool): Whether to skip files that already have outputs
    """
    
    # Ensure the output directory exists FIRST
    os.makedirs(output_summary_dir, exist_ok=True)
    
    # Set up logging AFTER creating directory
    logger = setup_logging(output_summary_dir)
    
    # Validate inputs
    if not os.path.exists(file_list_path):
        logger.error(f"File list not found: {file_list_path}")
        return
    
    # Read the file list CSV
    try:
        file_list_df = pd.read_csv(file_list_path)
        logger.info(f"Loaded file list from {file_list_path}: {len(file_list_df)} files to process")
    except Exception as e:
        logger.error(f"Error reading file list: {e}")
        return
    
    # Validate required columns
    required_columns = ['output_name', 'input_dir']
    missing_columns = [col for col in required_columns if col not in file_list_df.columns]
    if missing_columns:
        logger.error(f"Missing required columns in file list: {missing_columns}")
        return

    # Filter out already processed files if requested
    if skip_existing:
        to_process_df, already_processed_df = filter_processed_files(file_list_df, output_summary_dir)
        
        if len(already_processed_df) > 0:
            logger.info(f"Skipping {len(already_processed_df)} files that already have outputs:")
            for _, row in already_processed_df.iterrows():
                logger.info(f"  - {row['output_name']} (outputs exist)")
        
        if len(to_process_df) == 0:
            logger.info("All files have already been processed. Nothing to do!")
            return
            
        file_rows = [row for _, row in to_process_df.iterrows()]
        logger.info(f"Processing {len(file_rows)} remaining files")
    else:
        file_rows = [row for _, row in file_list_df.iterrows()]
        logger.info(f"Processing all {len(file_rows)} files (not checking for existing outputs)")
    
    logger.info(f"Starting parallel processing with max_workers={max_workers}")
    
    # Track results
    successful = 0
    failed = 0
    
    # Use ProcessPoolExecutor to parallelize file processing
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Submit all jobs
        future_to_file = {
            executor.submit(process_file, row, output_summary_dir, config_path, verbose): row['output_name']
            for row in file_rows
        }
        
        # Process results as they complete
        for future in as_completed(future_to_file):
            output_name = future_to_file[future]
            try:
                result = future.result()
                if result.startswith("SUCCESS"):
                    successful += 1
                    logger.info(f"✓ {result}")
                else:
                    failed += 1
                    logger.error(f"✗ {result}")
            except Exception as exc:
                failed += 1
                logger.error(f"✗ {output_name} generated an exception: {exc}")
    
    # Final summary
    total_in_list = len(file_list_df)
    skipped = len(already_processed_df) if skip_existing else 0
    processed = successful + failed
    
    logger.info(f"Processing complete!")
    logger.info(f"  Total files in list: {total_in_list}")
    if skip_existing:
        logger.info(f"  Skipped (already processed): {skipped}")
    logger.info(f"  Attempted to process: {processed}")
    logger.info(f"  Successful: {successful}")
    logger.info(f"  Failed: {failed}")
    
    if failed > 0:
        logger.warning(f"{failed} files failed to process. Check logs for details.")

def validate_file_list(file_list_path, limit_files=None):
    """
    Validate the file list and check if input directories exist.
    
    Parameters:
        file_list_path (str): Path to the file list CSV
        limit_files (int): Limit to first N files for testing (None for all)
        
    Returns:
        tuple: (valid_df, invalid_df) - DataFrames of valid and invalid entries
    """
    try:
        df = pd.read_csv(file_list_path)
        if limit_files is not None:
            df = df.head(limit_files)
            print(f"Limited to first {limit_files} files for testing")
    except Exception as e:
        print(f"Error reading file list: {e}")
        return None, None
    
    # Check if input directories exist  
    df['input_exists'] = df['input_dir'].apply(os.path.exists)
    
    valid_df = df[df['input_exists']].copy()
    invalid_df = df[~df['input_exists']].copy()
    
    print(f"File validation results:")
    print(f"  Total files: {len(df)}")
    print(f"  Valid files (input_dir exists): {len(valid_df)}")
    print(f"  Invalid files (input_dir missing): {len(invalid_df)}")
    
    if len(invalid_df) > 0:
        print(f"\nInvalid entries:")
        for _, row in invalid_df.iterrows():
            print(f"  {row['output_name']}: {row['input_dir']}")
    
    return valid_df, invalid_df

def check_existing_outputs(file_list_path, output_summary_dir, limit_files=None):
    """
    Check which files already have outputs and which need processing.
    
    Parameters:
        file_list_path (str): Path to the file list CSV
        output_summary_dir (str): Base output directory
        limit_files (int): Limit to first N files for testing (None for all)
    """
    try:
        df = pd.read_csv(file_list_path)
        if limit_files is not None:
            df = df.head(limit_files)
    except Exception as e:
        print(f"Error reading file list: {e}")
        return
    
    to_process, already_processed = filter_processed_files(df, output_summary_dir)
    
    print(f"Output status check:")
    print(f"  Total files: {len(df)}")
    print(f"  Already processed: {len(already_processed)}")
    print(f"  Need processing: {len(to_process)}")
    
    if len(already_processed) > 0:
        print(f"\nAlready processed:")
        for _, row in already_processed.iterrows():
            print(f"  ✓ {row['output_name']}")
    
    if len(to_process) > 0:
        print(f"\nNeed processing:")
        for _, row in to_process.iterrows():
            print(f"  - {row['output_name']}")

if __name__ == "__main__":
    # Configuration
    project_dir = "/mnt/data/malaria/synthetic_genomes/jessica_projects/FPG_ObsModelTesting"
    
    # Set the path for simulations from EMOD
    #file_list_path = os.path.join(project_dir, "6yr_itns_merged.csv")   
    file_list_path = os.path.join(project_dir, "test_mapping.csv")

    # Specify the directory for outputs
    output_summary_dir = os.path.join(project_dir, "infectionFPGReport_sympomaticsOnly")  
    
    # Optional config path (use None for defaults)
    config_path = None  # or specify path to your config.json
    
    # Check existing outputs first
    print("Checking existing outputs...")
    check_existing_outputs(file_list_path, output_summary_dir, limit_files=4)
    
    # Validate file list (optional but recommended)
    print("\nValidating file list...")
    valid_df, invalid_df = validate_file_list(file_list_path, limit_files=4)
    
    if valid_df is not None and len(valid_df) > 0:
        print(f"\nStarting processing of {len(valid_df)} valid files...")
        
        # Process files in parallel
        process_all_files(
            file_list_path=file_list_path,
            output_summary_dir=output_summary_dir,
            config_path=config_path,
            max_workers=6,  # Adjust based on your VM capabilities
            verbose=False,  # Set to True for detailed logging
            skip_existing=True  # Set to False to reprocess all files
        )
    else:
        print("No valid files to process!")
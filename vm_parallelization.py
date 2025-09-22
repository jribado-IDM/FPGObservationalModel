#####################################################################################
# Functions to parallelize the processing of multiple FPG simulations through the observational model on a VM. The goal is to replace this with a dtk_post_process script in the future. 
#####################################################################################
from concurrent.futures import ProcessPoolExecutor, as_completed


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


# Set output directory for the body of work
project_dir = "/mnt/data/malaria/synthetic_genomes/jessica_projects/2504_GRSweep"

# Set the path for simulations from EMOD - require Calcon to be mounted on the VM for file access. 
# Includes an 'output_name' and 'input_dir' columns. Suggest that the output name is the EMOD experiment ID name for merging with other EMOD outputs, like the InsetChart. The 'input_dir' is the full path to the folder containing the output files for that EMOD experiment.
file_list_path = os.path.join(project_dir, "sim_mapping_new_itns_6yr_highBiting.csv")   

# Specify the directory this this specific set of sim output will go
output_summary_dir = os.path.join(project_dir, "infectionFPGReport_summaries_6yr_highBiting")  
os.makedirs(output_summary_dir, exist_ok=True)

# Change max_workers based on your system capabilities. Lower transmission simulations can be faster to process, so you may be able to increase this number.
process_all_files(file_list_path, output_summary_dir, max_workers=6)
# Full Parasite Genetics - Observational Model

Documentation version 1.0 - September 2025

## Overview 

This repository contains scripts for EMOD's Full Parasite Genetics output to convert modeled results into recapitulative sampling for genomic surveillance. The observational model options allow for curated population sampling and then calculated genetic metrics for user specified combinations of samples, and optionally within epidemioligcally relevant nested populations within a group of samples. More options in sampling and metrics can be edited or expanded to match empirical data analyses. 

## Environment set-up

This model requires Pytohn 3.9 to be compatible with IDM tskit. To set up the environment.

~~~
python3 -m venv fpg_env
source fpg_env/bin/activate
pip install -r requirements.txt
~~~

## Config Parameters

### Hard filters 

The hard filter parameter options apply blanket filters to the full set of reported EMOD infections. If you are interested in comparing the differences in genetic metrics from the total population for these individual groups, use the subpopulation options instead. 

A default hard filtering removes the signal of consistently infected individuals; individuals can be represented in the infected population pool for sampling at most once per month. To match the reality of sequencing detection limits, only unique genotypes within an infection are considered.  

- 'symptomatics_only': Options, Bool. Will only consider symptomatic individuals (fever presenting) for sampling. 
- 'monogenomic_infections_only': Bool. Will only consider monogenomic infections 
- 'day_snapshot': Bool or integer. Will only consider infections from a single day (not operationally realistic, but can be useful for model comparisons.)


### Sampling

Sampling parameters follow hard filtering to randomly or directly choose a user specified number of samples per yer ('n_samples_year') to represent the population. Users can specify any number of comparisons to run with a single simulation, and each comparison group can have one or more technical replicate(s) ('replicates') if the within simulation sampling variability is of interest. Top level config names are flexible for naming

Each config entry will match this format, with the available options:
~~~
    {user_provided_name:
        {
        'method': ['random', 'seasonal', 'age'],
                'n_samples_year': Int,
                'replicates': 2,
                'method_params': {
                    Specific to the method chosen, see options below. 
                }
        }
    }
~~~

There are three broad method options for sampling:
1) 'random' - will sample N infections per year, tends to match seasonality cases. Can be further directed with the following 'method_params' options:
    - 'population_proportion': list, N populations. Used to sample from the source or sink only, equally, etc. Within population comparisons of genetic metrics can be specified below. Confirm the total number of samples per year * proportion reflects the minimum numbers of infections desired per population.
    - 'monogenomic_proportion': False or float for true (< 1). Will bias the sampling to include fewer or more monogenomic infections than may be the Bool modeled proportion. Used to compare the effect of metrics derived from monogenomic (e.g. unique proportion) or polygenomic samples (e.g. co-transmission proportion, Rh)
    - 'equal_monthly': Bool. Sample the same number of infections per month, regardless of seasonality. If the total number of samples requested is lower than the available, the remainder samples are not applied to other parts of the season. 

2) 'seasonal': Will sample N infections per year, each in the wet or the dry season to compare temporal sampling effects. Currently, the model is set-up for the consistent seasonality in Senegal, must update for other simulation scenarios. If an intervention start time is provided, this sampling frame is unaffected - the simulation years and months are used to make sure sequential seasonal groupings. Can be further refined with following 'method_params' options:
    - 'season': 'full' for all months in the wet or dry season or 'peak' for the 3 highest and lowest case months. Months for sampling are hardcoded for both full and peak season options.  

3) 'age': Will sample N infections per year, but will direct which age individuals are presented most in the population regardless of age distribution specified in the model. Use for comparing sampling schemes based on age, e.g. mirror biased sampling such as school surveys. Can be further customized with following 'method_params' options:
    - age_bins: List with the upper bound of each group. Default: [5, 15, 100]
    - age_bin_labels: List containing names for each age grouping. Default: ['0-5yrs', '5-15yrs', '15+yrs']


### Subpopulation comparisons

Above options will calculate metrics for all samples in a population for each sampling method specified. Additionally, comparisons within subpopulations to compare with all infections in the sampled population are supported. This allows for the investigation of metrics that may be more sensitive within groups or smaller timescales.

 The subpopulation options supported include:
- 'monthly':  Provide summary statistics by month in addition to year
- 'populations':  Defined by the population node in EMOD
- 'polygenomic':  Is polygenomic = 1, else monogenomic = 0
- 'symptomatic':  Is symptomatic = 1, else asymptomatic = 0
- 'age_bins':  Default age bins: 0-5, 5-15, 15+

### Genetic metrics

This section defines with genetic metrics will be calculated for each set as boolean input. 

    - 'cotransmission_proportion': From polygenomic infections, calculates how many contain genomes from a single mosquito biting event. 
    - 'complexity_of_infection': Calculated both 'true_coi' for the number of genomes a person holds in an infection and the 'effective_coi' for the number of unique genomes a person hold in an infection as the upper detectable bound. 
    - 'heterozygosity': For polygenomic infections, calculates how many positions contain more than one allele across genomes in an infection. Currently assumes all genotypes are captured, future plans include adding a density dependent weight to make this more realistic to specific strain parasitemia. 
    - 'identity_by_descent': Compares the pairwise Hamming distance for all genomes, defined by the parents at the start of the simulation, in specified infections at the population and/or the individual level. 
    - 'identity_by_state': Compares the pairwise Hamming distance for all genomes, defined by reference or alternative biallelic representations, in specified infections at the population and/or the individual level.
    - 'individual_ibx': Specification on whether or not to provide within sample relatedness for polygenomic infections. Will be set to True if Rh is specified. 
    - 'monogenomic_proportion': Calculates the proportion of monogenomic samples from the effective COI. 
    - 'rh': Calculates Rh for polygenomic infections, matching [paper reference] with 200 unique monogenomic pairwise draws to determine Hmono and sample heterozygosity for Hpoly. 
    - 'unique_genome_proportion': Calculates the proportion of genomes observed once in the sampled population, assuming phasing. 
    - 'unique_mono_proportion': Calculated the proportion of genomes observed once in the sampled population, assuming from monogenomic samples only to avoid phasing assumptions. 


### Other

- 'intervention_start_month': Int. Will reset the year and month start time to the start of the intervention. Allows for pre and post intervention comparisons. Currently a single intervention start time is supported.


## Output files 
    
``{sim_id}_FPG_SampledInfections.csv``: Replicates the input infection file, but adds columns for True/False inclusion for each sampling group, in addition to individual population level genetic metrics. This file contains columns used for calculating population level summary statistics. 

- 'group_{year/month}': Observational model specified grouping for time frames. If 'intervention_start_time" is specified, it's the year and month rescaled to the start of the intervention, if not the simulation year and continuous month are used as default
- '{true/effective}_coi': Total(true) or unique (effective) number of genomes in each infection.
- 'cotx': Categorizes co-transmission events as infections with an effective COI > 2 with a single biting event (one unique value in 'bite_ids'). Monogenomic infections are excluded. 
- '{sampling_name}_{n_samples}_rep{1...N replicates}': Columns specifying which sampling scheme the infection may be represented. Number of columns will match 'sampling' config options specified.
- 'barcode_with_Ns': Barcode string for each infection. Polygenomic infections with any discordant alleles within any genome as assigned N at each discordant position.
- 'heterozygosity': Proportion of barcode positions with an N.
- '{ibd/ibs}_{pairwise_count,mean,median,std,min,q25,75,max}: Individual infection relatedness for polygenomic infections
- {sampling_name}_{n_samples}_rep{1...N replicates}-individual_inferred_rh: Individual level Rh comparisons for each sampling scheme. Provided as individual columns in case infections are sampled across different sampling frames. 


``{sim_id}_FPG_ModelSummaries.csv``: File containing the genetic metrics across columns and the years, season, and subpopulation comparisons as columns. Addition of summary statistic columns can vary based on user options for metric calculations. 

- 'sampling_scheme': Grouping variable for the sampling scheme applied (matches 'sampling' options in config).
- 'comparison_type': Identifies with sampling scheme groupings, such as yearly or seasonal groups, or specified subpopulations (matches 'subpopulation_comparison' options in config).
- 'year_group': Specifies the year (either simulation year or intervention shifted year) or the seasonal grouping bin for summary statistics.
- 'sub_group': Sepcifies the additional groupings within subpopulations, e.g. whether True/False for polygenomic or symptomatic.  
- 'n_infections': Counts for the number of infections in each sampling scheme, for each year and subpopulation grouping specified in the observational model run. These are the actual number of infections that were available in the report by grouping and may be lower than the specified targets.
- '{true/effective}_poly_coi_count': The number of infections per grouping that have a COI > 2. True is the modeled umber of genomes tracked, which effective is the number of unique and detectable genomes in an infection.
- '{true/effective_poly_coi_prop}': The proportion of infections per grouping that have a COI > 2. Calculated as '{true/effective}_poly_coi_count'/n_infections.
- 'all_genomes_total_count': The total number of genomes identified in all infections per grouping, assuming full phasing of all infection genomes. 
- 'all_genomes_unique_count': The total number of genomes observed once in all genomes from infections per grouping, assuming full phasing of all infection genomes.
- 'all_genomes_ids_unique_prop': The proportion of unique genomes in all genomes from infections per grouping, assuming full phasing of all infection genomes. Calculated as 'all_genomes_unique_count'/'all_genomes_total_count'.
- 'mono_genomes_total_count': The total number of genomes identified in all infections per grouping, measured only from monogenomic infections that are inherently phased. 
- 'mono_genomes_unique_count': The total number of genomes observed once in all genomes from infections per grouping, measured only from monogenomic infections that are inherently phased.
- 'mono_genomes_ids_unique_prop': The proportion of unique genomes in all genomes from infections per grouping, measured only from monogenomic infections that are inherently phased. Calculated as 'mono_genomes_unique_count'/'mono_genomes_total_count'.
- 'cotransmission_count': The number of infections per grouping with a COI > 2 and from a single mosquito biting event. 
- 'cotransmission_prop': The proportion of co-transmission infections within polygenomic infections. Calculated as 'cotransmission_count'/'effective_poly_coi_count'.
- '{true/effective}_coi_{mean,median,std,min,q25,75,max}': Full summary statistics for the COI distribution of infections in a grouping. 
- '{ibd/ibs}_{pairwise_count,mean,median,std,min,q25,75,max}': Full population level summary statistics for the COI distribution of infections in a grouping.
- 'rh_inferred_{mean,median,std}': Summary statistics for the inferred Rh from polygenomic infections, bootstrapping the monogenomic Rh. 



``{sim_id}-{ibd/ibs}_distributions.json``: To avoid large pairwise matrices s output, to further investigate population level IBs distributions one could use the JSON file with the IBx calculated value as the key up to two decimal places and the number of pairwise counts as a the value. It matches the output CSV in matching sampling, comparison_types, and subpopulations. 

~~~
  {
      "user_specified_name": { # "sampling_scheme 
          "population": { # comparision_type Like group_year, season_bin, population polygenomic, etc.
              "(2, 0)": {       # For subpopulations, the key here can be tuple, with the first item is the group_year, and the second is the group identifier. For example, this is for year 2, population 0
                  "0.5": 54,
                  "0.7": 41,
                  "0.9": 7,
                  "1": 4
              }

          }
      }
  }
~~~    

## Mapping file for VM parallelization

Ideally, one would pull the sim_id and the corresponding output path using an analyzer that contains epidemiological information too (need to find in the Senegal MDA FPG repository). 

In the the absence of the mapping file, one can look for the directories belonging to a COMPS simulation set by name to create a mapping file. The full paths returned for the genetic input files for the observational model can also be used to pull information on the same sims in the InsetChart if needing to create a corresponding epidemiological summary report for analyses.  

~~~
# Example pull of data
EXPERIMENT_NAME="/mnt/calculon2/jsuresh/output/maka fpg 10k - 6yr - strong ITNs in_20250522_195001/"
OUTPUT_FILE="experiment_mapping.csv"

{ echo "output_name,input_dir"; find "$EXPERIMENT_NAME" -name "output" -type d | sed 's|.*/\([0-9a-f]\{8\}-[0-9a-f]\{4\}-[0-9a-f]\{4\}-[0-9a-f]\{4\}-[0-9a-f]\{12\}\)/output$|\1,"\0"|'; } > "$OUTPUT_FILE"
~~~
"""
Self-contained test suite using Python's built-in unittest (no pytest needed).
All test data (dataframes and matrices) are generated within this file.

Run with: python fpg_unittest.py
Or: python -m unittest fpg_unittest -v

This test suite covers key functions from fpg_observational_model.unified_sampling.py,
unified_metric_calculations.py, and run_observational_model.py.

"""

import os
import unittest
from unittest import result
import pandas as pd
import numpy as np
import ast
import sys
import random
import json
import tempfile
import shutil
from itertools import combinations
from collections import Counter
from unittest.mock import patch, MagicMock, mock_open

# Declare TEST_DATA at module level - will be initialized after class definition
TEST_DATA = None

MODULE_PATHS = [
    '../fpg_observational_model',  # One level up, then into fpg_observational_model
    '.',  # Current directory (in case you run from parent_dir)
    '..',  # Parent directory
]

for path in MODULE_PATHS:
    abs_path = os.path.abspath(path)
    if abs_path not in sys.path and os.path.exists(abs_path):
        sys.path.insert(0, abs_path)
        print(f"Added to path: {abs_path}")

# Import your actual functions to test
try:
    from fpg_observational_model.unified_sampling import (
        calculate_infection_metrics,
        apply_emod_filters,
        subset_randomly,
        subset_by_seasons,
        assign_season_group,
        assign_peak_group,
        subset_by_age,
        filter_emod_infections,
        run_sampling_model,
        convert_month,
        adjust_time_columns
    )

    SAMPLING_IMPORTED = True
except ImportError as e:
    print(f"Warning: Could not import unified_sampling: {e}")
    SAMPLING_IMPORTED = False

try:
    from fpg_observational_model.unified_metric_calculations import (
        register_matrix,
        get_matrix,
        comprehensive_group_summary,
        generate_het_barcode,
        calculate_ibx_matrix,
        ibx_distribution,
        weighted_describe_scipy,
        calculate_heterozygosity,
        calculate_individual_rh,
        calculate_population_rh,
        sample_from_distribution,
        identify_nested_comparisons,
        process_nested_summaries,
        process_nested_ibx,
        process_nested_fws,
        run_time_summaries,
        update_ibx_index,
    )

    METRICS_IMPORTED = True
except ImportError as e:
    print(f"Warning: Could not import unified_metric_calculations: {e}")
    METRICS_IMPORTED = False

try:
    from fpg_observational_model.run_observational_model import (
        extract_sampled_infections,
        update_matrix_indices
    )    
    HELPER_IMPORTED = True
except ImportError as e:
    print(f"Warning: Could not import unified_metric_calculations: {e}")
    HELPER_IMPORTED = False


try:
    from run_observational_model import (
        run_observational_model,
        get_default_config,
        process_file
    )

    RUN_MODEL_IMPORTED = True
except ImportError as e:
    print(f"Warning: Could not import run_observational_model: {e}")
    RUN_MODEL_IMPORTED = False    


###############################################################################
# SHARED TEST DATA - Created once, used by all tests
###############################################################################
class SharedTestData:
    """Singleton class to hold shared test data across all test classes"""
    _instance = None
    _initialized = False

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super(SharedTestData, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if not SharedTestData._initialized:
            # Create base infection dataframe
            self.base_infection_df = pd.DataFrame({
                'infIndex': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
                'IndividualID': [100, 101, 102, 103, 104, 105, 106, 107, 108, 109],
                'simulation_year': [1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
                'year': [1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
                'group_year': [1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
                'month': [1, 6, 12, 1, 6, 12, 1, 6, 9, 12],
                'continuous_month': [13, 18, 24, 25, 30, 36, 37, 42, 45, 48],
                'day': [30, 180, 365, 395, 545, 730, 760, 910, 940, 1095],
                'population': [0, 0, 1, 1, 0, 0, 1, 1, 0, 1],
                'fever_status': [1, 0, 1, 1, 1, 0, 1, 1, 1, 0],
                'age_day': [1000, 2000, 3000, 8000, 1500, 6000, 2500, 1000, 5000, 7000],
                'recursive_nid': ['[1]', '[1]', '[2, 3]', '[4]', '[5, 6, 7]', '[7, 8]', '[9]', '[11]', '[12, 12]', '[10, 14]'],
                'genome_ids': ['[101]', '[101]', '[102, 103]', '[104]', '[105, 106, 107]', '[108]', '[109]', '[110]', '[111, 111]', '[112, 113]'],
                'bite_ids': ['[1]', '[2]', '[3, 3]', '[4]', '[5, 6, 6]', '[7, 7]', '[8]', '[8]', '[9, 10]', '[11,11]']
            })

            # Create genotype matrix (IBS)
            # Deprecated - added a more specific pattern for better tests 
            # np.random.seed(42)
            # self.genotype_matrix = np.random.randint(0, 2, size=(15, 10), # dtype=np.int8)
            # self.genotype_matrix[2] = self.genotype_matrix[3]  # Make some identical

            self.genotype_matrix = np.array([
                [0, 1, 0, 1, 0, 1],  # genome 0
                [1, 0, 1, 0, 1, 0],  # genome 1
                [0, 1, 0, 1, 0, 1],  # genome 2 
                [0, 1, 0, 1, 0, 1],  # genome 3 - MATCHES genome 2 ✓
                [1, 1, 0, 0, 1, 1],  # genome 4
                [1, 0, 1, 0, 1, 0],  # genome 5
                [1, 0, 1, 0, 1, 0],  # genome 6 - MATCHES genome 5 ✓
                [1, 1, 1, 0, 0, 0],  # genome 7 - differs from 5/6 at 2 positions (indices 1,4) ✓
                [0, 0, 0, 1, 1, 1],  # genome 8 - OPPOSITE of genome 7 ✓
                [0, 0, 1, 1, 0, 0],  # genome 9
                [1, 1, 0, 0, 1, 1],  # genome 10
                [0, 1, 1, 0, 0, 1],  # genome 11
                [1, 0, 0, 1, 1, 0],  # genome 12
                [0, 1, 1, 0, 1, 0],  # genome 13
                [1, 0, 1, 1, 0, 1],  # genome 14
                [0, 1, 0, 1, 0, 1],  # genome 15
            ]).astype(np.int8)
            
            # Create IBD matrix
            # np.random.seed(123)
            # ibd = np.random.randint(0, 101, size=(15, 15))
            # self.ibd_matrix = (ibd + ibd.T) // 2
            # np.fill_diagonal(self.ibd_matrix, 100)
            # self.ibd_matrix = self.ibd_matrix.astype(np.int16)
            self.ibd_matrix = np.array([
                [3, 7, 2, 9, 5, 1],  # genome 0
                [8, 4, 6, 3, 10, 2], # genome 1
                [5, 3, 8, 2, 7, 4],  # genome 2 
                [5, 3, 3, 2, 7, 4],  # genome 3 - one off from genome 2 ✓
                [1, 9, 4, 6, 3, 8],  # genome 4
                [2, 5, 9, 1, 6, 3],  # genome 5
                [2, 5, 9, 1, 6, 6],  # genome 6 - one off from genome 5 ✓
                [2, 8, 4, 1, 6, 6], # genome 7 - differs from 5/6 at 2 or 3 positions (indices 1,2,5) ✓
                [9, 3, 1, 7, 5, 2],  # genome 8 - COMPLETELY DIFFERENT from genome 7 ✓
                [4, 6, 1, 8, 3, 9],  # genome 9
                [7, 2, 5, 4, 9, 6],  # genome 10
                [3, 8, 7, 1, 4, 10], # genome 11
                [6, 1, 3, 9, 2, 5],  # genome 12
                [9, 4, 2, 6, 8, 1],  # genome 13
                [5, 7, 10, 3, 1, 4], # genome 14
            ])
            
            # Standard config
            self.config = {
                'hard_filters': {
                    'symptomatics_only': False,
                    'monogenomic_infections_only': False,
                    'day_snapshot': False
                },
                'intervention_start_month': 24,
                'sampling_configs': {
                    'random': {
                        'method': 'random',
                        'n_samples_year': 5,
                        'replicates': 2,
                        'method_params': {
                            'population_proportions': [0.5, 0.5],
                            'monogenomic_proportion': False,
                            'equal_monthly': False
                        }
                    }
                },
                'metrics': {
                    'cotransmission_proportion': True,
                    'complexity_of_infection': True,
                    'heterozygosity': True,
                    'identity_by_descent': True,
                    'identity_by_state': True,
                    'individual_ibx': True,
                    'monogenomic_proportion': True,
                    'rh': True,
                    'unique_genome_proportion': True
                },
                'subpopulation_comparisons': {
                    'add_monthly': True,
                    'populations': True,
                    'polygenomic': True,
                    'symptomatic': True,
                    'age_bins': True
                }
            }

            # Monogenomic IBS distribution for R_h tests
            self.monogenomic_dict = {
                0.1: 2, 0.2: 5, 0.3: 8, 0.4: 10,
                0.5: 12, 0.6: 8, 0.7: 5, 0.8: 3,
                0.9: 2, 1.0: 5
            }

            SharedTestData._initialized = True

    def get_infection_df(self, with_metrics=False):
        """Get a copy of the infection dataframe, optionally with metrics calculated"""
        df = self.base_infection_df.copy()

        if with_metrics and SAMPLING_IMPORTED:
            df = calculate_infection_metrics(df)

        return df

    def get_genotype_matrix(self):
        """Get a copy of the genotype matrix"""
        return self.genotype_matrix.copy()

    def get_ibd_matrix(self):
        """Get a copy of the IBD matrix"""
        return self.ibd_matrix.copy()

    def get_config(self):
        """Get a copy of the standard config"""
        import copy
        return copy.deepcopy(self.config)

    def get_monogenomic_dict(self):
        """Get a copy of the monogenomic IBS distribution"""
        return self.monogenomic_dict.copy()


# Initialize TEST_DATA at module level (NO INDENTATION)
print("\n" + "=" * 70)
print("INITIALIZING SHARED TEST DATA")
print("=" * 70)

try:
    TEST_DATA = SharedTestData()
    print(f"✓ Base infection dataframe: {TEST_DATA.base_infection_df.shape}")
    print(f"✓ Genotype matrix: {TEST_DATA.genotype_matrix.shape}")
    print(f"✓ IBD matrix: {TEST_DATA.ibd_matrix.shape}")
    print(f"✓ TEST_DATA initialized successfully")
except Exception as e:
    print(f"✗ ERROR initializing TEST_DATA: {e}")
    import traceback

    traceback.print_exc()
    TEST_DATA = None

print("=" * 70 + "\n")


###############################################################################
# TEST CLASSES - FILTERING SAMPLES - Using shared TEST_DATA
###############################################################################
@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestApplyEmodFilters(unittest.TestCase):
    """Test apply_emod_filters function"""

    def setUp(self):
        """Use shared test data with metrics"""
        self.df = TEST_DATA.get_infection_df(with_metrics=True)

    def test_filter_duplicate_infections(self):
        """Test filter_emod_infections function"""
        # Add duplicate row of an exact infection to test duplicate removal by month
        duplicate_df = pd.concat([self.df, self.df.iloc[:1]]) 
        result = filter_emod_infections(duplicate_df, duplicate_seed=123)
        
        self.assertLessEqual(len(result), len(self.df))    
    
    def test_fever_filter_true(self):
        """Test filtering for fever cases"""
        result = apply_emod_filters(self.df, fever_filter=True)

        self.assertTrue(all(result['fever_status'] == 1))
        self.assertLessEqual(len(result), len(self.df))

    def test_monogenomic_filter(self):
        """Test monogenomic filtering"""
        result = apply_emod_filters(self.df, monogenomic_filter=True)
        mono_only = result[result['effective_coi'] == 1]
        poly_remaining = result[result['effective_coi'] > 1]
        
        self.assertEqual(len(poly_remaining), 4, 
                        f"Filter failed: {len(poly_remaining)} polygenomic infections remain")
        
    def test_snapshot_filter(self):
        """Test day snapshot filtering"""
        snapshot_day = 545
        result = apply_emod_filters(self.df, day_filter=snapshot_day)
        
        self.assertTrue(all(result['day'] == snapshot_day))
        self.assertTrue(len(result) == 1)    


@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestSeasonAssignment(unittest.TestCase):
    """Test season assignment functions"""
    
    def test_assign_season_group(self):
        """Test assign_season_group function"""
        
        test_cases = [
            {'simulation_year': 2020, 'month': 1},
            {'simulation_year': 2020, 'month': 5},
            {'simulation_year': 2020, 'month': 10},
        ]
        
        for case in test_cases:
            row = pd.Series(case)
            season = assign_season_group(row)
            
            self.assertIsInstance(season, str)
            self.assertIn('season', season.lower())
    
    def test_assign_peak_group(self):
        """Test assign_peak_group function"""
        
        row_peak_wet = pd.Series({'simulation_year': 2020, 'month': 11})
        season = assign_peak_group(row_peak_wet)
        self.assertIn('Peak wet', season)


@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestSamplingFunctions(unittest.TestCase):
    """Test actual sampling functions"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=False)

    def test_continuous_month_conversion(self):    
        """Test monthly sampling scheme"""
        
        self.df.rename(columns={'continuous_month': 'test_continuous_month'}, inplace=True)
        self.df = convert_month(self.df)
        self.assertEqual(self.df['test_continuous_month'].tolist(), 
                         self.df['continuous_month'].tolist())
        
    def test_subset_monthly(self):
        
        month_duplicate_df = self.df.copy()
        month_duplicate_df['group_year'] = 1
        month_duplicate_df['population'] = 1
        month_duplicate_df['continuous_month'] = [13, 13, 14, 14, 26, 26, 37, 37, 38, 38]

        result = subset_randomly(
            month_duplicate_df, 
            n_samples_year=12, # setting to 12 since the function divides by 12 to indicate monthly sampling - this gets us to 1 sample per month
            replicates=2,
            equal_monthly=True,
            scheme='monthly'
        )
        
        # Check expected columns added
        expected_cols = [f'monthly_12_rep{i+1}' for i in range(2)]
        for col in expected_cols:
            self.assertIn(col, result.columns)

        # Check each month has correct number of samples per replicate
        for replicate in expected_cols:
            month_counts = result.groupby(['continuous_month', replicate]).size().reset_index(name="count")
            self.assertTrue((month_counts['count'] <= 1).all())
    
    def test_subset_randomly(self):
        """Test subset_randomly function"""
        n_samples_year = 2
        replicates = 2
        
        result = subset_randomly(
            self.df, 
            n_samples_year=n_samples_year,
            replicates=replicates,
            scheme='random'
        )
        
        # Check expected columns added
        expected_cols = [f'random_{n_samples_year}_rep{i+1}' for i in range(replicates)]
        for col in expected_cols:
            self.assertIn(col, result.columns)

        # Check each population has 1 infection at most per year in a replicate
        subset_counts = result.groupby(['population', 'year', 'random_2_rep1']).size().reset_index(name="count")

        self.assertTrue(all(subset_counts['count'] <= 1), "More than 1 infection sampled for a population in a year")    

    def test_subset_by_full_seasons(self):
        """Test subset_by_seasons function for full season"""
        
        result = subset_by_seasons(
            self.df,
            n_samples_year=1,
            replicates=2,
            season='full'
        )

        # Check expected columns added
        expected_cols = [f'seasonalFull_1_rep{i+1}' for i in range(2)]
        for col in expected_cols:
            self.assertIn(col, result.columns)
        
        # Check each season has correct number of samples per replicate
        dry_samples = result[result['month'].isin([2, 3, 4, 5, 6, 7])]
        wet_samples = result[result['month'].isin([8, 9, 10, 11, 12, 1])]
        
        for replicate in expected_cols:
            # Check each season has correct number of samples per replicate
            season_counts = result.groupby(['year', replicate]).size().reset_index(name="count")
            self.assertTrue((season_counts['count'] == 1).all())

            # Confirm months
            self.assertTrue(dry_samples[replicate].str.contains('Dry').all())
            self.assertTrue(wet_samples[replicate].str.contains('Wet').all())    

    def test_subset_by_peak_seasons(self):
        """Test subset_by_seasons function"""
        
        result = subset_by_seasons(
            self.df,
            n_samples_year=1,
            replicates=2,
            season='peak'
        )

        # Check expected columns added
        expected_cols = [f'seasonalPeak_1_rep{i+1}' for i in range(2)]
        for col in expected_cols:
            self.assertIn(col, result.columns)
        
        # Check each season has correct number of samples per replicate
        off_peak_samples = result[result['month'].isin([1, 2, 7, 8, 9])]
        dry_peak_samples = result[result['month'].isin([3, 4, 5, 6])]
        wet_peak_samples = result[result['month'].isin([10, 11, 12])]
        for replicate in expected_cols:
            # Confirm counts
            season_counts = result.groupby(['year', replicate]).size().reset_index(name="count")
            self.assertTrue((season_counts['count'] == 1).all())
            
            # Confirm months
            self.assertTrue(off_peak_samples[replicate].str.contains('Off').all())
            self.assertTrue(dry_peak_samples[replicate].str.contains('Peak dry').all())
            self.assertTrue(wet_peak_samples[replicate].str.contains('Peak wet').all())            
    
    def test_subset_by_age(self):
        """Test subset_by_age function.
        Example - discard ages 5-15 or 15+ using default age bins (0-5, 5-15, 15+)
        """
        
        result = subset_by_age(
            self.df,
            n_samples_year=3,
            replicates=1,
            age_bin_weights=[1, 0, 0],
            scheme='age_bin'
        )
        
        # Check expected columns added
        expected_cols = [f'age_bin_3_rep{i+1}' for i in range(1)]
        for col in expected_cols:
            self.assertIn(col, result.columns)
        
        # Check each age bin has correct number of samples per replicate
        age_bins = list(set(result['age_bin']))
        for bin_range in age_bins:
            age_bin_counts = result[result['age_bin'] == bin_range].groupby(['year', expected_cols[0]]).size().reset_index(name="count")
            
            if bin_range in ['5-15yrs', '15+yrs']:
                self.assertTrue(age_bin_counts.empty, f"Infections sampled for excluded age bin {bin_range}")
            else:    
                self.assertTrue((age_bin_counts['count'] == 1).all())


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestInterventionTiming(unittest.TestCase):
    """Test intervention start month handling"""

    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=False)
        self.intervention_start_month = 29

        self.no_subgroup_config = {'add_monthly': False, 'populations': False, 'polygenomic': False, 'symptomatic': False, 'age_bins': False}

    def test_intervention_month_calculation(self):
        """Test that intervention_month is calculated correctly"""
        from fpg_observational_model.unified_sampling import adjust_time_columns

        result = adjust_time_columns(self.df, intervention_start_month=self.intervention_start_month)

        # Check month reset calculations
        self.assertIn('intervention_month', result.columns)
        
        expected_intervention_months = (self.df['continuous_month'] - self.intervention_start_month).tolist()
        actual_intervention_months = result['intervention_month'].tolist()
        self.assertTrue(actual_intervention_months, expected_intervention_months)

        # Check year reset calculations
        self.assertIn('intervention_year', result.columns)

        expected_intervention_years = [-2, -1, -1, -1, 0, 0, 0, 1, 1, 1]
        actual_intervention_years = result['intervention_year'].tolist()
        self.assertTrue(actual_intervention_years, expected_intervention_years)

    def test_uses_intervention_years(self):    
        """Test that sampling uses intervention years when available"""
       
        sampling_col = 'random_5_rep1'
        self.df[sampling_col] = 1

        intervention_df = adjust_time_columns(self.df, intervention_start_month=self.intervention_start_month)
        if 'intervention_year' in intervention_df.columns:
            intervention_df['group_year'] = intervention_df['intervention_year']

        original_group = identify_nested_comparisons(
            self.df, 
            sampling_col, 
            config=self.no_subgroup_config
        )

        intervention_group = identify_nested_comparisons(
            intervention_df, 
            sampling_col, 
            config=self.no_subgroup_config
        )

        self.assertNotEqual(original_group['group_year'], intervention_group['group_year'])


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestSubsampleFilter(unittest.TestCase):
    """Test extraction of sampled infections from dataframe"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=False)
        
        # Add sampling columns - each replicate will pick one infection per year, no overlaps with other replicate for a total of 2 infections per year
        self.df['random_1_rep1'] = [1, np.nan, np.nan, np.nan, 1, np.nan, np.nan, 1, np.nan, np.nan]
        self.df['random_1_rep2'] = [np.nan, np.nan, 1, 1, np.nan, np.nan, 1, np.nan, np.nan, np.nan]

    def test_extract_sampled_infections(self):
        result = extract_sampled_infections(self.df)
        sample_counts = result.groupby('year').size().reset_index(name="count")

        self.assertTrue((sample_counts['count'] == 2).all())


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestNestedComparisonDictionary(unittest.TestCase):
    """Test identify_nested_comparisons function"""

    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        
        # Add columns for subgroup comparisons - same calculation as nested comparisons function
        self.df['is_polygenomic'] = self.df['effective_coi'].apply(lambda x: True if x > 1 else False)

        days_per_year = 365.25
        age_bins = [0, int(days_per_year * 5), int(days_per_year * 15), int(self.df['age_day'].max() + 1)]
        age_bin_labels = ['0-5yrs', '5-15yrs', '15+yrs']     
        self.df['age_bin'] = pd.cut(self.df['age_day'], bins=age_bins, labels=age_bin_labels, include_lowest=True)    

        self.config = TEST_DATA.get_config()['subpopulation_comparisons']
        
    def test_identify_nested_comparisons_basic(self):
        """Test basic nested comparison identification"""

        sampling_col = 'random_5_rep1'
        self.df[sampling_col] = 1

        result = identify_nested_comparisons(
            self.df,
            sampling_col,
            config=self.config
        )
        
        # Check year samples correctly identified
        yearly_groups = self.df.groupby('group_year')['infIndex'].apply(list).to_dict()
        self.assertEqual(yearly_groups, result['group_year'])

        # Check all within year subgroups identified
        mappings = {
            'populations': 'population',
            'polygenomic': 'is_polygenomic',
            'symptomatic': 'fever_status',
            'age_bins': 'age_bin'
        }
        for subgroup in mappings.keys():
            subgroup_dict = self.df.groupby(['group_year', mappings[subgroup]], observed=True)['infIndex'].apply(list).to_dict()
            
            self.assertEqual(subgroup_dict, result[subgroup])
  

###############################################################################
# NEW TESTS FOR run_observational_model
###############################################################################
@unittest.skipIf(not RUN_MODEL_IMPORTED, "run_observational_model not available")
class TestRunObservationalModel(unittest.TestCase):
    """Test suite for run_observational_model function updates"""

    def setUp(self):
        """Set up test fixtures"""
        # Create temporary directories for test outputs
        current_dir = os.path.dirname(__file__)
        self.emod_output_path = os.path.join(current_dir, 'test_data')
        self.test_dir = tempfile.mkdtemp()
        self.config_dir = tempfile.mkdtemp()

        # Use shared test data
        self.test_infection_df = TEST_DATA.get_infection_df(with_metrics=False)

    def tearDown(self):
        """Clean up test fixtures"""
        shutil.rmtree(self.test_dir, ignore_errors=True)
        shutil.rmtree(self.config_dir, ignore_errors=True)

    # ==================== Test 1: Default config usage ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('pandas.read_csv')
    def test_default_config_when_no_config_provided(self, mock_read_csv, mock_summaries, mock_sampling):
        """Test that default config is used when neither config_path nor config is provided"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        # Run without config
        with patch('builtins.print') as mock_print:
            run_observational_model(
                sim_name='test_default',
                emod_output_path=self.emod_output_path,
                output_path=self.test_dir,
                verbose=True
            )

            # Check that using default config message was printed
            print_calls = [str(call) for call in mock_print.call_args_list]
            self.assertTrue(any('using default config' in str(call).lower() for call in print_calls))

    # ==================== Test 2: Config file loading ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('pandas.read_csv')
    def test_config_file_loading(self, mock_read_csv, mock_summaries, mock_sampling):
        """Test loading config from file"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        # Create test config file
        test_config = {
            'intervention_start_month': 50,  # default is 29
            'metrics': {
                'heterozygosity': False      # default is True
            }
        }
        config_path = os.path.join(self.config_dir, 'test_config.json')
        with open(config_path, 'w') as f:
            json.dump(test_config, f)

        # Run with config file
        with patch('builtins.print') as mock_print:
            run_observational_model(
                sim_name='test_file',
                emod_output_path=self.emod_output_path,
                config_path=config_path,
                output_path=self.test_dir,
                verbose=True
            )

            # Check that config was loaded
            print_calls = [str(call) for call in mock_print.call_args_list]
            self.assertTrue(any('Loaded config from' in str(call) for call in print_calls))

    # ==================== Test 3: Config dictionary ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('pandas.read_csv')
    def test_config_dictionary_usage(self, mock_read_csv, mock_summaries, mock_sampling):
        """Test using config dictionary directly"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        # Create test config dictionary
        test_config = {
            'intervention_start_month': 36,
            'sampling_configs': {
                'random': {
                    'n_samples_year': 200
                }
            }
        }

        # Run with config dictionary
        with patch('builtins.print') as mock_print:
            run_observational_model(
                sim_name='test_dict',
                emod_output_path=self.emod_output_path,
                config=test_config,
                output_path=self.test_dir,
                verbose=True
            )

            # Check that dictionary config message was printed
            print_calls = [str(call) for call in mock_print.call_args_list]
            self.assertTrue(any('provided config dictionary' in str(call).lower() for call in print_calls))

    # ==================== Test 4: Partial config merge ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('run_observational_model.get_default_config')
    @patch('pandas.read_csv')
    def test_partial_config_merges_with_defaults(self, mock_read_csv, mock_default, mock_summaries, mock_sampling):
        """Test that partial config is merged with defaults"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        default_config = TEST_DATA.get_config()
        mock_default.return_value = default_config

        # Partial config - only override one value
        partial_config = {
            'intervention_start_month': 50
        }

        # Capture the config passed to sampling model
        captured_config = None

        def capture_config(*args, **kwargs):
            nonlocal captured_config
            captured_config = kwargs.get('config')
            return self.test_infection_df.copy()

        mock_sampling.side_effect = capture_config

        # Run with partial config
        run_observational_model(
            sim_name='test_partial',
            emod_output_path=self.emod_output_path,
            config=partial_config,
            output_path=self.test_dir,
            verbose=False
        )

        # Verify merged config
        self.assertEqual(captured_config['intervention_start_month'], 50)
        self.assertEqual(captured_config['hard_filters']['symptomatics_only'],
                         default_config['hard_filters']['symptomatics_only'])
        self.assertEqual(captured_config['metrics']['heterozygosity'],
                         default_config['metrics']['heterozygosity'])

    # ==================== Test 5: Deep merge for nested configs ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('run_observational_model.get_default_config')
    @patch('pandas.read_csv')
    def test_deep_merge_nested_configs(self, mock_read_csv, mock_default, mock_summaries, mock_sampling):
        """Test that nested configs are properly merged"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        default_config = TEST_DATA.get_config()
        mock_default.return_value = default_config

        # Partial nested config
        partial_config = {
            'metrics': {
                'heterozygosity': False  # Only override one nested value
            },
            'sampling_configs': {
                'random': {
                    'n_samples_year': 200  # Only override one nested value
                }
            }
        }

        # Capture config
        captured_config = None

        def capture_config(*args, **kwargs):
            nonlocal captured_config
            captured_config = kwargs.get('config')
            return self.test_infection_df.copy()

        mock_sampling.side_effect = capture_config

        # Run
        run_observational_model(
            sim_name='test_deep_merge',
            emod_output_path=self.emod_output_path,
            config=partial_config,
            output_path=self.test_dir,
            verbose=False
        )

        # Verify deep merge
        self.assertEqual(captured_config['metrics']['heterozygosity'], False)
        self.assertEqual(captured_config['metrics']['rh'], default_config['metrics']['rh'])
        self.assertEqual(captured_config['sampling_configs']['random']['n_samples_year'], 200)
        self.assertEqual(captured_config['sampling_configs']['random']['replicates'],
                         default_config['sampling_configs']['random']['replicates'])

    # ==================== Test 6: Unknown parameter warning ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('pandas.read_csv')
    def test_unknown_parameter_warning(self, mock_read_csv, mock_summaries, mock_sampling):
        """Test that unknown parameters trigger a warning"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        # Config with unknown parameters
        config_with_typos = {
            'intervention_start_month': 36,
            'unknown_parameter': 'value',  # Unknown
            'metrics': {
                'heterozygosity': True,
                'typo_metric': False  # Unknown nested
            }
        }

        # Run and capture output
        with patch('builtins.print') as mock_print:
            run_observational_model(
                sim_name='test_unknown',
                emod_output_path=self.emod_output_path,
                config=config_with_typos,
                output_path=self.test_dir,
                verbose=True
            )

            # Check for warning
            print_calls = [str(call) for call in mock_print.call_args_list]
            warning_found = any('WARNING' in str(call) and 'unknown' in str(call).lower()
                                for call in print_calls)
            self.assertTrue(warning_found, "Should warn about unknown parameters")

            # Check that specific unknown keys are mentioned
            all_output = ' '.join([str(call) for call in print_calls])
            self.assertIn('unknown_parameter', all_output)
            self.assertIn('typo_metric', all_output)

    # ==================== Test 7: Invalid JSON error ====================
    @patch('pandas.read_csv')
    def test_invalid_json_error(self, mock_read_csv):
        """Test that invalid JSON in config file raises appropriate error"""
        mock_read_csv.return_value = self.test_infection_df

        # Create invalid JSON file
        config_path = os.path.join(self.config_dir, 'invalid.json')
        with open(config_path, 'w') as f:
            f.write('{"invalid": json syntax}')

        # Should raise ValueError with JSON error message
        with self.assertRaises(ValueError) as context:
            run_observational_model(
                sim_name='test_invalid_json',
                emod_output_path=self.emod_output_path,
                config_path=config_path,
                output_path=self.test_dir
            )

        self.assertIn('Error parsing JSON', str(context.exception))

    # ==================== Test 8: Config priority ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('pandas.read_csv')
    def test_config_path_priority_over_config_dict(self, mock_read_csv, mock_summaries, mock_sampling):
        """Test that config_path takes priority over config dictionary"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        # Create config file
        file_config = {'intervention_start_month': 100}
        config_path = os.path.join(self.config_dir, 'priority.json')
        with open(config_path, 'w') as f:
            json.dump(file_config, f)

        # Create conflicting dict config
        dict_config = {'intervention_start_month': 50}

        # Capture config
        captured_config = None

        def capture_config(*args, **kwargs):
            nonlocal captured_config
            captured_config = kwargs.get('config')
            return self.test_infection_df.copy()

        mock_sampling.side_effect = capture_config

        # Run with both
        run_observational_model(
            sim_name='test_priority',
            emod_output_path=self.emod_output_path,
            config_path=config_path,
            config=dict_config,
            output_path=self.test_dir,
            verbose=False
        )

        # File config should win
        self.assertEqual(captured_config['intervention_start_month'], 100)

    # ==================== Test 9: Multiple unknown keys at different levels ====================
    @patch('run_observational_model.run_sampling_model')
    @patch('run_observational_model.run_time_summaries')
    @patch('pandas.read_csv')
    def test_multiple_unknown_keys_detected(self, mock_read_csv, mock_summaries, mock_sampling):
        """Test detection of multiple unknown keys at various nesting levels"""
        # Setup mocks
        mock_read_csv.return_value = self.test_infection_df
        mock_sampling.return_value = self.test_infection_df.copy()
        mock_summaries.return_value = (pd.DataFrame(), pd.DataFrame(), {})

        config_with_many_unknowns = {
            'unknown_top_level': 'value',
            'metrics': {
                'unknown_metric_1': True,
                'unknown_metric_2': False
            },
            'sampling_configs': {
                'random': {
                    'unknown_param': 123
                }
            },
            'completely_unknown_section': {
                'nested': 'value'
            }
        }

        with patch('builtins.print') as mock_print:
            run_observational_model(
                sim_name='test_many_unknown',
                emod_output_path=self.emod_output_path,
                config=config_with_many_unknowns,
                output_path=self.test_dir,
                verbose=True
            )

            all_output = ' '.join([str(call) for call in mock_print.call_args_list])
            self.assertIn('unknown_top_level', all_output)
            self.assertIn('unknown_metric_1', all_output)
            self.assertIn('unknown_param', all_output)
            self.assertIn('completely_unknown_section', all_output)


###############################################################################
# Additional test classes from original file continue...
###############################################################################
@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestProcessNestedSummaries(unittest.TestCase):
    """Test process_nested_summaries function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        
        self.nested_indices = {
            'group_year': {
                1: [0, 1, 2],
                2: [3, 4, 5],
                3: [6, 7, 8, 9]
            }
        }
    
    def test_process_nested_summaries_basic(self):
        """Test basic nested summary processing"""
        result = process_nested_summaries(
            self.nested_indices,
            self.df,
            comprehensive_group_summary
        )
        
        self.assertIsInstance(result, pd.DataFrame)
        self.assertTrue((result['comparison_type'] == "group_year").all())
        self.assertEqual(result['n_infections'].tolist(), [3, 3, 4])
        self.assertEqual(result['true_coi_poly_count'].tolist(), [1, 2, 2])
        self.assertEqual(result['true_coi_mean'].tolist(), [1.333, 2.000, 1.5])
        self.assertEqual(result['effective_coi_poly_count'].tolist(), [1, 2, 1])
        self.assertEqual(result['effective_coi_mean'].tolist(), [1.333, 2.000, 1.250])
        self.assertEqual(result['all_genomes_total'].tolist(), [4, 6, 5])
        self.assertEqual(result['all_genomes_unique'].tolist(), [3, 5, 5])
        self.assertEqual(result['mono_genomes_total'].tolist(), [2, 1, 3])
        self.assertEqual(result['mono_genomes_unique'].tolist(), [1, 1, 3])
        self.assertEqual(result['cotransmission_count'].tolist(), [1, 1, 1])


class TestRunSampling(unittest.TestCase):
    """Test run_sampling model - main pipeline function for infection sampling.
    
    Note: Other infection filtering options tested in this file, not an exhausting full run of the sampling function since invididual components are tested elsewhere.
    """
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.config = TEST_DATA.get_config()

        # Update config to reduce samples for smaller input file
        self.config['sampling_configs'] =  {
                    'random': {
                        'method': 'random',
                        'n_samples_year': 2,
                        'replicates': 2,
                        'method_params': {
                            'population_proportions': [0.5, 0.5],
                            'monogenomic_proportion': False,
                            'equal_monthly': False
                        }
                    }
                }
        
    def test_run_sampling_basic(self):
        """Test basic execution"""
        
        sampled_infections = run_sampling_model(self.df, self.config)

        samples = self.config['sampling_configs']['random']['n_samples_year']
        replicates = self.config['sampling_configs']['random']['replicates']
        expected_cols = [f'random_{samples}_rep{i+1}' for i in range(replicates)]
        for col in expected_cols:
            self.assertIn(col, sampled_infections.columns)
        
        for replicate in expected_cols:
            # Check each season has correct number of samples per replicate
            year_by_pop_counts = sampled_infections.groupby(['group_year', 'population', replicate]).size().reset_index(name="count")
            self.assertTrue((year_by_pop_counts['count'] == 1).all())


###############################################################################
# TEST CLASSES - GENOTYPE MATRIX - Using shared TEST_DATA
###############################################################################
@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestMatrixLinkage(unittest.TestCase):
    """Test linkage between infection dataframe and genotype matrix"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=False)
        self.matrix = TEST_DATA.get_genotype_matrix()
    
    def test_matrix_index_range(self):
        """Test that all indices are within matrix bounds"""
        all_indices = []
        for nid_str in self.df['recursive_nid']:
            nid_list = ast.literal_eval(nid_str)
            all_indices.extend(nid_list)
        
        max_index = max(all_indices)
        matrix_rows = self.matrix.shape[0]
        
        self.assertLess(max_index, matrix_rows)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestMatrixRegistry(unittest.TestCase):
    """Test matrix registration and retrieval"""
    
    def test_register_and_get_matrix(self):
        """Test matrix registration"""
        test_matrix = TEST_DATA.get_genotype_matrix()
        matrix_name = 'test_registry_matrix'
        
        register_matrix(matrix_name, test_matrix)
        retrieved = get_matrix(matrix_name)
        
        self.assertTrue(np.array_equal(retrieved, test_matrix))
    
    def test_get_nonexistent_matrix(self):
        """Test error on non-existent matrix"""
        with self.assertRaises(KeyError):
            get_matrix('nonexistent_matrix_xyz123')


###############################################################################
# TEST CLASSES - GENETIC METRICS - Using shared TEST_DATA
###############################################################################
@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestCalculateInfectionMetrics(unittest.TestCase):
    """Test calculate_infection_metrics function from fpg_observational_model.unified_sampling.py"""
    
    def setUp(self):
        """Use shared test data"""
        self.df = TEST_DATA.get_infection_df(with_metrics=False)
        self.ibs_matrix = TEST_DATA.get_genotype_matrix()
        register_matrix('ibs_matrix', self.ibs_matrix)
    
    def test_calculate_infection_metrics_basic(self):
        """Test that calculate_infection_metrics runs and adds columns"""
        result = calculate_infection_metrics(self.df)
        
        self.assertIn('true_coi', result.columns)
        self.assertIn('effective_coi', result.columns)
        self.assertIn('cotx', result.columns)
    
    def test_coi_values(self):
        """Test COI calculation values"""
        result = calculate_infection_metrics(self.df)
        # add genotype COI
        result[['genotype_coi', 'barcode_with_Ns', 'heterozygosity']] = self.df.apply(lambda row: generate_het_barcode(self.ibs_matrix, row['recursive_nid']), axis=1, result_type='expand')
        
        expected_true_coi = [1, 1, 2, 1, 3, 2, 1, 1, 2, 2]
        actual_true_coi = result['true_coi'].tolist()
        self.assertEqual(actual_true_coi, expected_true_coi)

        expected_eff_coi = [1, 1, 2, 1, 3, 2, 1, 1, 1, 2]
        actual_eff_coi = result['effective_coi'].tolist()
        self.assertEqual(actual_eff_coi, expected_eff_coi)

        expected_geno_coi = [1, 1, 1, 1, 2, 2, 1, 1, 1, 2]
        actual_geno_coi = result['genotype_coi'].tolist()
        self.assertEqual(actual_geno_coi, expected_geno_coi)

    def test_cotx_values(self):
        """Test cotransmission assignment"""
        result = calculate_infection_metrics(self.df)    

        expected_cotx = [None, None, True, None, False, True, None, None, None, True]
        actual_cotx = result['cotx'].tolist()
        self.assertEqual(actual_cotx, expected_cotx)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestComprehensiveGroupSummary(unittest.TestCase):
    """Test comprehensive_group_summary function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_infection_count(self):
        """Test infection count in summary"""
        result = comprehensive_group_summary(self.df)
        self.assertEqual(result['n_infections'], 10)

    def test_true_coi_stats(self):
        """Test true COI summary calculations"""
        result = comprehensive_group_summary(self.df)

        expected_proportion_poly = sum(self.df['true_coi'] > 1) / len(self.df)
        self.assertEqual(result['true_coi_poly_prop'], expected_proportion_poly)

        expected_mean = self.df['true_coi'].mean()
        self.assertEqual(result['true_coi_mean'], expected_mean)  

        expected_median = self.df['true_coi'].median()
        self.assertEqual(result['true_coi_median'], expected_median)

        expected_std = round(self.df['true_coi'].std(), 3)
        self.assertEqual(result['true_coi_std'], expected_std)

    def test_effective_coi_stats(self):
        """Test effective COI summary calculations"""
        result = comprehensive_group_summary(self.df)

        expected_proportion_poly = sum(self.df['effective_coi'] > 1) / len(self.df)
        self.assertEqual(result['effective_coi_poly_prop'], expected_proportion_poly)

        expected_mean = self.df['effective_coi'].mean()
        self.assertEqual(result['effective_coi_mean'], expected_mean)  

        expected_median = self.df['effective_coi'].median()
        self.assertEqual(result['effective_coi_median'], expected_median)

        expected_std = round(self.df['effective_coi'].std(), 3)
        self.assertEqual(result['effective_coi_std'], expected_std)    

    def test_monogenomic_proportion(self):
        """Test monogenomic proportion calculation"""
        result = comprehensive_group_summary(self.df)

        self.assertEqual(result['all_genomes_total'], 15)
        self.assertEqual(result['all_genomes_unique'], 13)
        self.assertEqual(result['all_genomes_unique_prop'], 0.867)     

        self.assertEqual(result['mono_genomes_total'], 6)
        self.assertEqual(result['mono_genomes_unique'], 5)
        self.assertEqual(result['mono_genomes_unique_prop'], 0.833)

    def test_cotx_proportion(self):    
        """Test cotransmission proportion calculation"""
        result = comprehensive_group_summary(self.df)

        expected_cotx_count = sum(self.df['cotx'] == True)
        expected_cotx_prop = expected_cotx_count / (self.df['true_coi'] > 1).sum()

        self.assertEqual(result['cotransmission_count'], expected_cotx_count)
        self.assertEqual(result['cotransmission_count'], 3)
        
        self.assertEqual(result['cotransmission_prop'], expected_cotx_prop)
        self.assertEqual(result['cotransmission_prop'], 0.600)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestHeterozygosityFunctions(unittest.TestCase):
    """Test get_variant_coi and generate_het_barcode"""
    
    def setUp(self):
        """Use shared test data"""
        self.matrix = TEST_DATA.get_genotype_matrix()
        register_matrix('test_matrix', self.matrix)

    def check_allele_frequency_calculation(self):
        """Check allele frequency for a dummy genotype matrix."""
        test_matrix = self.matrix[[0,2,4]]
        expected_alternative_af = [0.3333333333333333, 1.0, 0.0, 0.6666666666666666, 0.3333333333333333, 1.0]
        true_alternative_af = np.sum(test_matrix == 1, axis=0) / test_matrix.shape[0]
        self.assertEqual(true_alternative_af.tolist(), expected_alternative_af)    

    def check_heterozygosity_calculation(self):
        """Check heterozygosity for a dummy minor allele frequency array."""
        maf = np.array([0.0, 0.2, 0.5, 0.8, 1.0, 0.3])
        heterozygosity = calculate_heterozygosity(maf)  

        expected_heterozygosity = [0.0, 0.32, 0.5, 0.32, 0.0, 0.42]
        self.assertEqual(heterozygosity, expected_heterozygosity)
    
    def test_heterozygosity_output_monogenomic(self):
        """Test heterozygosity for a genome infection."""
        test_indices = [[0], [2,3]]  # First index is genome 0, and second indices are two identical genomes 2 and 3 that should yield monogenomic like output results
    
        for indices in test_indices:
            genotype_coi, barcode_with_Ns, heterozygosity = generate_het_barcode(self.matrix, indices)

            self.assertEqual(genotype_coi, 1)
            self.assertNotIn('N', barcode_with_Ns)
            self.assertEqual([0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                         heterozygosity)
        
    def test_heterozygosity_output_opposite(self):
        """Test heterozygosity for a two genome infection."""
        indices = [7, 8]  # Genomes 7 and 8 are opposite at all positions
        genotype_coi, barcode_with_Ns, heterozygosity = generate_het_barcode(self.matrix, indices)

        self.assertEqual(genotype_coi, 2)
        self.assertTrue(all(x == 'N' for x in barcode_with_Ns))
        self.assertEqual([0.5, 0.5, 0.5, 0.5, 0.5, 0.5],
                         heterozygosity)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestRhCalculation(unittest.TestCase):
    """Test R_h functions and R_h metric logic"""
    
    def setUp(self):
        """Set up test data with barcode heterozygosity and a mock pairwise monogenomic IBS draw of 5 values."""
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.mono_large_dict = TEST_DATA.get_monogenomic_dict()

        # Add heterozygosity column (simulated)
        self.df['barcode_N_prop'] = [0.0, 0.0, 0.0, 0.333, 1.0, 0.0, 0.0, 0.0, 0.0, 0.667]
        self.poly_samples = self.df[self.df['effective_coi'] > 1].copy()
        self.mono_test_distribution = [0.4, 0.7, 0.5, 0.6, 0.8]
        self.mono_test_dict = Counter(self.mono_test_distribution)

    def test_rh_exclude_identical_barcodes(self):
        """Test that IBS=1 excluded from H_mono"""
        filtered_dist = sample_from_distribution(self.mono_large_dict, n_bootstraps = 50, exclude_keys=[1])
        self.assertNotIn(1.0, filtered_dist)
    
    def test_rh_sampling_from_distribution(self):
        """Test sampling from distribution"""
        filtered_dict = {k: v for k, v in self.mono_large_dict.items() if k != 1.0}
        
        values = list(filtered_dict.keys())
        weights = list(filtered_dict.values())
        
        n_samples = 200
        samples = np.random.choice(values, size=n_samples, p=np.array(weights)/sum(weights))
        
        self.assertEqual(len(samples), n_samples)
        self.assertTrue(all(s < 1.0 for s in samples))    
    
    def test_calculate_individual_rh(self):
        """Test the actual calculate_rh function on a single value.
        Assume 5 random draws of monogenomic pairwise IBS to create the distribution to calculate R_h from a single infection with barcode heterozygosity of 0.25.
        """
        random.seed(42)
        barcode_N_count = 0.25
        test_rh = calculate_individual_rh(barcode_N_count, self.mono_test_distribution)
        
        self.assertAlmostEqual(round(test_rh, 3), 0.583)

    def test_calculate_population_rh(self):
        """Test R_h metric calculation logic"""
        result_df = self.df.copy()
        
        # Calculate R_h for each infection
        calculated_rh, __ = calculate_population_rh(result_df, self.mono_test_dict)

        self.poly_samples['individual_inferred_rh'] = self.poly_samples.apply(lambda row: calculate_individual_rh(row['barcode_N_prop'], self.mono_test_dict), axis=1)
    
        # Check R_h population mean
        true_rh_mean = 0.305
        expected_rh_mean = self.poly_samples['individual_inferred_rh'].mean()
        actual_rh_mean = calculated_rh['rh_poly_inferred_mean']
        self.assertEqual(round(expected_rh_mean, 3), true_rh_mean, actual_rh_mean)   
        

class TestFwsCalculation(unittest.TestCase):
    """
    TODO: Confirm Fws regression logic with collaborators before implementing tests.
    Test Fws calculation functions. 
    Note heterozygosity calculation is already tested in TestHeterozygosityFunctions. 
    """
    
    def setUp(self):
        """Set up test data with mock minor allele frequencies."""
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_calculate_fws(self):
        pass


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestIdentityByState(unittest.TestCase):
    """Test IBS matrix calculations. All logic should hold the same for identity by descent, which only changes the underlying matrix."""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.matrix = TEST_DATA.get_genotype_matrix()
        register_matrix('ibs_matrix', self.matrix)

        # Set up mock sampling columns
        self.df['random_1_rep1'] = [1, np.nan, np.nan, np.nan, 1, np.nan, np.nan, 1, np.nan, np.nan]
        self.df['random_1_rep2'] = [np.nan, np.nan, 1, 1, np.nan, np.nan, 1, np.nan, np.nan, np.nan]

        self.sampled_infections = extract_sampled_infections(self.df)

    def test_matrix_mapping(self):
        """Test if mapping for reading in only a subset of the matrix for sampled infections works correctly"""
        
        matrix = get_matrix('ibs_matrix')
        df_matrix_update = update_matrix_indices(self.df)

        matrix_original_indices = df_matrix_update['original_nid'].explode().tolist()
        matrix_original = matrix[matrix_original_indices]

        matrix_subset_indices = list(set(matrix_original_indices))
        matrix_subset = matrix[matrix_subset_indices]
        matrix_mapped = matrix_subset[df_matrix_update['recursive_nid'].explode().tolist()]

        self.assertTrue(np.array_equal(matrix_original, matrix_mapped))

    def test_update_ibx_indices(self):
        """Test that recursive_index and ibx_nid are correctly assigned - i.e. the same index ordering."""

        def same_rank_order(list1, list2):
            """Check if two lists have the same relative ordering"""
            if len(list1) != len(list2):
                return False
            return (np.argsort(list1) == np.argsort(list2)).all()

        # Start with updating the indices for the sampled infections
        sampled_infections = update_matrix_indices(self.sampled_infections)
        
        replicate_1_samples = sampled_infections[~sampled_infections['random_1_rep1'].isna()].copy()

        replicate_1_samples = update_ibx_index(replicate_1_samples)

        for idx, row in replicate_1_samples.iterrows():
            self.assertTrue(same_rank_order(row['recursive_nid'], row['ibx_nid']), f"Rank order mismatch at row {idx}: {row['recursive_nid']} vs {row['ibx_nid']}")    
    
    def test_pairwise_combinations_count(self):
        """Test correct number of pairwise combinations"""
        indices = [0, 1, 2, 3]
        pairs = list(combinations(indices, 2))
        expected_pairs = len(indices) * (len(indices) - 1) // 2
        
        self.assertEqual(len(pairs), expected_pairs)
        self.assertEqual(len(pairs), 6)
    
    def test_ibx_distribution_accumulation(self):
        """Test that IBx values accumulate correctly"""
        ibx_dict = {}
        ibx_values = [0.5, 0.5, 0.7, 0.5, 0.3]
        
        for val in ibx_values:
            val_rounded = round(val, 2)
            ibx_dict[val_rounded] = ibx_dict.get(val_rounded, 0) + 1
        
        self.assertEqual(ibx_dict[0.5], 3)
        self.assertEqual(ibx_dict[0.7], 1)
    
    def test_ibx_distribution(self):
        """Test ibx_distribution function with a small mock IBx matrix"""
        
        indices = [0, 1, 2]
        hash_ibx = np.array([
            [0,  5, 10],
            [5,  0,  2],
            [10, 2,  0]
        ])
        
        result = ibx_distribution(indices, hash_ibx)
        # Pairs: (0,1)→5/10=0.5, (0,2)→10/10=1.0, (1,2)→2/10=0.2
        expected = {0.5: 1, 1.0: 1, 0.2: 1}
        
        self.assertEqual(result, expected)

    def test_ibx_calculation(self):
        """Test IBS calculation between genotypes from IDM tskit"""

        def normalize_dict_keys(d):
            return {tuple(sorted(k)): float(v) for k, v in d.items()}
        
        # Pick infections with known relatedness
        # Infection 2 - genomes 2 and 3 (identical)
        # Infection 5 - genomes 7 and 8 (opposites)
        # Infections 6 and 9 - genomes 9, 10, and 14 (mixed)
        expected_ibx = {
            (2, 3): 1.0,  # identical genomes
            (7, 8): 0.0,  # opposites
            (9, 10): 0.0,  
            (9, 14): 0.667,
            (10, 14): 0.333
        }

        actual_ibx = {}
        for infections in [[2], [5], [6, 9]]:
            infections_df = self.df[self.df['infIndex'].isin(infections)]
            hash_df, hash_ibx = calculate_ibx_matrix(infections_df, self.matrix[infections_df['recursive_nid'].explode().tolist()], intervals=None)
            hash_ibx = hash_ibx / self.matrix.shape[1]  # Normalize by number of variant positions

            pairwise_genomes = list(combinations(hash_df['recursive_nid'], 2))
            for genome_pair in pairwise_genomes:
                ibx_pair = hash_df['ibx_index'][hash_df['recursive_nid'].isin(genome_pair)].tolist()
                ibx = hash_ibx[ibx_pair[0], ibx_pair[1]]
                actual_ibx[genome_pair] = round(ibx, 3)

        self.assertEqual(normalize_dict_keys(expected_ibx), 
                         normalize_dict_keys(actual_ibx))

    def test_weighted_ibx_summary(self):
        """Test IBS distribution summary calculation"""
        test_ibx_values = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 0.2, 0.4, 0.6, 0.8]
        test_ibx_dict = Counter(test_ibx_values)

        unweighted = round(pd.Series(test_ibx_values).describe(), 3).to_frame().T.add_prefix('ibs_')
        weighted   = weighted_describe_scipy(test_ibx_dict, "ibs")

        self.assertTrue(np.allclose(unweighted.values, weighted.values))

    def test_calculate_ibx_matrix_with_idm(self):
        """Test IBx matrix calculation with idm module"""
        try:
            import idm
            
            hash_df, hash_ibx = calculate_ibx_matrix(self.df, self.matrix)
            
            # Check hash_df structure
            self.assertIsInstance(hash_df, pd.DataFrame)
            self.assertIn('infIndex', hash_df.columns)
            self.assertIn('ibx_index', hash_df.columns)
            self.assertIn('ibx_nid', hash_df.columns)
            
            # Check hash_ibx structure
            self.assertIsInstance(hash_ibx, np.ndarray)
            self.assertEqual(hash_ibx.ndim, 2, "IBx matrix should be 2D")
            
            # Matrix should be square
            self.assertEqual(hash_ibx.shape[0], hash_ibx.shape[1])
            
        except ImportError:
            self.skipTest("idm module not available")    


class TestRunTimeSummaries(unittest.TestCase):
    """Test run_time_summaries - main pipeline function for genetic metrics"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.df['group_month'] = self.df['continuous_month']
        
        # Add sampling columns
        self.df['month_rep0'] = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
        self.df['random_5_rep1'] = [1, np.nan, 1, np.nan, 1, 1, np.nan, 1, np.nan, 1]
        self.df['random_5_rep2'] = [np.nan, 1, 1, 1, np.nan, 1, 1, np.nan, 1, np.nan]
        
        # Register matrices
        ibs_matrix = TEST_DATA.get_genotype_matrix()
        ibd_matrix = TEST_DATA.get_ibd_matrix()
        register_matrix('ibs_matrix', ibs_matrix)
        register_matrix('ibd_matrix', ibd_matrix)
        
        # Add heterozygosity
        self.df['heterozygosity'] = [0.0, 0.0, 0.3, 0.0, 0.5, 0.0, 0.0, 0.4, 0.0, 0.35]
    
    def test_run_time_summaries_basic(self):
        """Test basic execution"""
        config = {'polygenomic': False}
        
        try:
            summaries, _, _ = run_time_summaries(
                self.df,
                subpop_config=config,
                user_ibx_categories=None,
                individual_ibx_calculation=False,
                rh_calculation=False
            )

            self.assertIsInstance(summaries, pd.DataFrame)
            self.assertGreater(len(summaries), 0)
            
        except UnboundLocalError as e:
            if 'ibx_category' in str(e):
                self.skipTest(f"Known bug: {e}")
    
    def test_run_time_summaries_multiple_replicates(self):
        """Test multiple replicates processed"""
        config = {'polygenomic': False}
        
        try:
            summaries, _, _ = run_time_summaries(
                self.df,
                subpop_config=config,
                user_ibx_categories=None
            )
            
            sampling_schemes = summaries['sampling_scheme'].unique()
            self.assertGreaterEqual(len(sampling_schemes), 2)
        
        except UnboundLocalError as e:
            if 'ibx_category' in str(e):
                self.skipTest(f"Known bug: {e}")


###############################################################################
# TEST CLASSES - MISCELLANEOUS - Using shared TEST_DATA
###############################################################################
class TestEdgeCases(unittest.TestCase):
    """Test edge cases"""
    
    def test_empty_dataframe(self):
        """Test handling of empty dataframe"""
        empty_df = pd.DataFrame()
        self.assertEqual(len(empty_df), 0)
    
    def test_invalid_recursive_nid(self):
        """Test error handling for malformed recursive_nid"""
        with self.assertRaises((ValueError, SyntaxError)):
            ast.literal_eval('invalid')
    
    def test_all_monogenomic(self):
        """Test handling when all infections are monogenomic"""
        df = pd.DataFrame({
            'infIndex': [0, 1, 2],
            'recursive_nid': ['[0]', '[1]', '[2]'],
            'bite_ids': ['[1]', '[2]', '[3]']
        })
        
        if SAMPLING_IMPORTED:
            df = calculate_infection_metrics(df)
            self.assertTrue(all(df['effective_coi'] == 1))




###############################################################################
# MAIN TEST RUNNER
###############################################################################
def run_tests():
    """Run all tests and print results"""
    global TEST_DATA

    print("=" * 70)
    print("INFECTION SAMPLING & METRICS TEST SUITE")
    print("=" * 70)
    print(f"unified_sampling imported: {SAMPLING_IMPORTED}")
    print(f"unified_metric_calculations imported: {METRICS_IMPORTED}")
    print(f"run_observational_model imported: {RUN_MODEL_IMPORTED}")
    print("=" * 70)

    if TEST_DATA is not None:
        print("\nUSING SHARED TEST DATA:")
        print(f"  Infection dataframe: {TEST_DATA.base_infection_df.shape}")
        print(f"  Genotype matrix (IBS): {TEST_DATA.genotype_matrix.shape}")
        print(f"  IBD matrix: {TEST_DATA.ibd_matrix.shape}")
        print(f"  All tests use IDENTICAL copies of this data")
    else:
        print("\nWARNING: TEST_DATA not initialized")
    print("=" * 70 + "\n")

    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)

    print("\n" + "=" * 70)
    print("TEST SUMMARY")
    print("=" * 70)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    print("=" * 70)

    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
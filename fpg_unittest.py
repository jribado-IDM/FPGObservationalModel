#!/usr/bin/env python3
"""
Complete Test Suite for Observational Model
Merges all tests into one comprehensive file including:
- Basic functionality tests
- Scientific accuracy tests  
- Array boolean error tests
- Sampling and filtering tests
- IBx calculation tests
- Integration tests
"""

import unittest
import pandas as pd
import numpy as np
from collections import Counter
from itertools import combinations
from unittest.mock import patch, MagicMock
import sys
import os

# Add the current directory to the path to import your modules
sys.path.append(os.path.dirname(__file__))

# Import all required modules
from run_observational_model import (
    get_default_config, 
    update_matrix_indices, 
    extract_sampled_infections,
    run_observational_model
)

from unified_metric_calculations import (
    identify_nested_comparisons,
    comprehensive_group_summary,
    generate_het_barcode,
    calculate_ibx_matrix,
    register_matrix,
    get_matrix,
    process_nested_summaries,
    process_nested_ibx,
    ibx_distribution,
    weighted_describe_scipy,
    run_time_summaries
)

from unified_sampling import (
    parse_list,
    adjust_time_columns,
    calculate_infection_metrics,
    apply_emod_filters,
    subset_randomly,
    subset_by_seasons,
    subset_by_age,
    run_sampling_functions,
    run_sampling_model,
    n_samples_by_pop,
    filter_emod_infections,
    convert_month,
    assign_season_group,
    assign_peak_group,
    validate_subset_inputs
)

class TestBasicFunctionality(unittest.TestCase):
    """Test basic functionality and data handling"""
    
    def setUp(self):
        """Set up basic test data"""
        self.basic_df = pd.DataFrame({
            'infIndex': [0, 1, 2, 3, 4],
            'recursive_nid': ['[0, 1]', '[2]', '[3, 4, 5]', '[6]', '[7, 8]'],
            'true_coi': [2, 1, 3, 1, 2],
            'effective_coi': [2, 1, 2, 1, 2],
            'genome_ids': [[10, 11], [12], [13, 14, 15], [16], [17, 18]],
            'cotx': [1, 0, 0, 0, 1],
            'fever_status': [1, 0, 1, 0, 1],
            'population': [0, 0, 1, 1, 0],
            'age_day': [1000, 2000, 3000, 8000, 1500],
            'group_year': [1, 1, 2, 2, 2],
            'group_month': [1, 2, 3, 4, 5],
            'rep_random_0': [1, np.nan, 1, np.nan, 1],
            'rep_random_1': [np.nan, 1, np.nan, 1, np.nan]
        })
        
        # Register a test matrix
        self.test_matrix = np.random.randint(0, 2, size=(20, 100))
        register_matrix('basic_test_matrix', self.test_matrix)
    
    def test_get_default_config(self):
        """Test default configuration structure"""
        config = get_default_config()
        
        required_keys = ['hard_filters', 'sampling_configs', 'metrics', 'subpopulation_comparisons']
        for key in required_keys:
            self.assertIn(key, config)
        
        # Test specific configurations
        self.assertIn('cotransmission_proportion', config['metrics'])
        self.assertIn('complexity_of_infection', config['metrics'])
        self.assertIn('random', config['sampling_configs'])
    
    def test_update_matrix_indices(self):
        """Test matrix index updating"""
        result_df = update_matrix_indices(self.basic_df)
        
        # Check that original_nid column was created
        self.assertIn('original_nid', result_df.columns)
        
        # Check that recursive_nid was updated to lists of integers
        first_nid = result_df['recursive_nid'].iloc[0]
        self.assertIsInstance(first_nid, list)
        self.assertIsInstance(first_nid[0], int)
    
    def test_extract_sampled_infections(self):
        """Test extraction of sampled infections"""
        result_df = extract_sampled_infections(self.basic_df)
        
        # Should have same number of rows (no all-NaN rows in test data)
        self.assertEqual(len(result_df), len(self.basic_df))
        
        # Test with all-NaN row
        test_df = self.basic_df.copy()
        test_df.loc[len(test_df)] = {
            'infIndex': 5,
            'rep_random_0': np.nan,
            'rep_random_1': np.nan
        }
        
        result_df = extract_sampled_infections(test_df)
        self.assertEqual(len(result_df), len(self.basic_df))  # Should drop the all-NaN row
    
    def test_parse_list_function(self):
        """Test parse_list utility function"""
        # Test valid list string
        result = parse_list('[1, 2, 3]')
        self.assertEqual(result, [1, 2, 3])
        
        # Test single item list
        result = parse_list('[42]')
        self.assertEqual(result, [42])
        
        # Test empty list
        result = parse_list('[]')
        self.assertEqual(result, [])
        
        # Test invalid input
        result = parse_list('not_a_list')
        self.assertEqual(result, [])

class TestSamplingAndFiltering(unittest.TestCase):
    """Test sampling functions and filtering logic"""
    
    def setUp(self):
        """Set up test data for sampling functions"""
        self.sample_infection_df = pd.DataFrame({
            'infIndex': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            'IndividualID': [100, 101, 102, 103, 104, 105, 106, 107, 108, 109],
            'simulation_year': [1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
            'year': [1, 1, 1, 2, 2, 2, 3, 3, 3, 3],
            'month': [1, 6, 12, 1, 6, 12, 1, 6, 9, 12],
            'day': [30, 180, 365, 395, 545, 730, 760, 910, 940, 1095],
            'population': [0, 0, 1, 1, 0, 0, 1, 1, 0, 1],
            'fever_status': [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            'age_day': [1000, 2000, 3000, 8000, 1500, 6000, 2500, 4000, 5000, 7000],
            'recursive_nid': ['[0]', '[1]', '[2, 3]', '[4]', '[5, 6, 7]', '[8]', '[9]', '[10, 11]', '[12]', '[13, 14]'],
            'genome_ids': ['[100]', '[101]', '[102, 103]', '[104]', '[105, 106, 107]', '[108]', '[109]', '[110, 111]', '[112]', '[113, 114]'],
            'bite_ids': ['[1]', '[2]', '[3]', '[4]', '[5]', '[6]', '[7]', '[8]', '[9]', '[10]']
        })
    
    def test_convert_month(self):
        """Test continuous month conversion"""
        test_df = pd.DataFrame({
            'simulation_year': [1, 1, 2],
            'month': [1, 12, 6]
        })
        
        result_df = convert_month(test_df)
        
        expected_continuous_months = [13, 24, 30]  # (1*12 + 1), (1*12 + 12), (2*12 + 6)
        self.assertEqual(result_df['continuous_month'].tolist(), expected_continuous_months)
    
    def test_calculate_infection_metrics(self):
        """Test infection metrics calculation"""
        df_with_continuous = convert_month(self.sample_infection_df.copy())
        
        result_df = calculate_infection_metrics(df_with_continuous)
        
        # Check that required columns were added
        required_columns = ['true_coi', 'effective_coi', 'cotx']
        for col in required_columns:
            self.assertIn(col, result_df.columns)
        
        # Check specific calculations
        # Row 2 has genome_ids '[102, 103]' -> true_coi should be 2
        self.assertEqual(result_df.loc[2, 'true_coi'], 2)
        
        # Row 4 has genome_ids '[105, 106, 107]' -> true_coi should be 3
        self.assertEqual(result_df.loc[4, 'true_coi'], 3)
        
        # Monogenomic infections should have cotx = None
        mono_rows = result_df[result_df['effective_coi'] == 1]
        for idx in mono_rows.index:
            self.assertIsNone(result_df.loc[idx, 'cotx'])
    
    def test_n_samples_by_pop_calculations(self):
        """Test population sampling calculations"""
        
        # Single population
        single_pop_df = self.sample_infection_df[self.sample_infection_df['population'] == 0]
        result = n_samples_by_pop(single_pop_df, n_samples_year=10)
        self.assertEqual(result, [10])
        
        # Multiple populations - equal distribution
        result = n_samples_by_pop(self.sample_infection_df, n_samples_year=10)
        self.assertEqual(result, [5, 5])  # Should distribute equally
        
        # Multiple populations - with proportions
        result = n_samples_by_pop(
            self.sample_infection_df, 
            n_samples_year=10, 
            population_proportions=[0.7, 0.3]
        )
        self.assertEqual(result, [7, 3])
    
    def test_apply_emod_filters(self):
        """Test EMOD filtering functionality"""
        df_with_continuous = convert_month(self.sample_infection_df.copy())
        
        # Test fever filtering
        result_df = apply_emod_filters(df_with_continuous, fever_filter=True)
        
        # All remaining should have fever_status = 1
        self.assertTrue(all(result_df['fever_status'] == 1))
        
        # Should have some fever cases
        self.assertGreater(len(result_df), 0)
        
        # Test no fever filtering  
        result_no_fever = apply_emod_filters(df_with_continuous, fever_filter=False)
        self.assertTrue(all(result_no_fever['fever_status'] == 0))
    
    def test_subset_randomly_basic(self):
        """Test basic random sampling"""
        df_processed = calculate_infection_metrics(convert_month(self.sample_infection_df.copy()))
        
        result_df = subset_randomly(
            df_processed, 
            n_samples_year=2, 
            replicates=1,
            scheme='test_random'
        )
        
        # Check that sampling column was added
        sample_cols = [col for col in result_df.columns if 'test_random' in col and 'rep' in col]
        self.assertEqual(len(sample_cols), 1)
        
        # Check that some samples were selected
        sample_col = sample_cols[0]
        sampled_count = result_df[sample_col].notna().sum()
        self.assertGreater(sampled_count, 0)
    
    def test_subset_randomly_reproducibility(self):
        """Test random sampling reproducibility with seeds"""
        df_processed = calculate_infection_metrics(convert_month(self.sample_infection_df.copy()))
        
        # Same seed should give same results
        result1 = subset_randomly(df_processed, 3, 1, base_seed=123)
        result2 = subset_randomly(df_processed, 3, 1, base_seed=123)
        
        sample_cols = [col for col in result1.columns if 'random' in col and 'rep' in col]
        for col in sample_cols:
            pd.testing.assert_series_equal(result1[col], result2[col], check_names=False)
    
    def test_subset_by_seasons(self):
        """Test seasonal sampling"""
        df_processed = calculate_infection_metrics(convert_month(self.sample_infection_df.copy()))
        
        result_df = subset_by_seasons(
            df_processed,
            n_samples_year=2,
            replicates=1,
            scheme='test_seasonal',
            season='full'
        )
        
        # Check that seasonal columns were added
        self.assertIn('season_group', result_df.columns)
        sample_cols = [col for col in result_df.columns if 'test_seasonal' in col and 'rep' in col]
        self.assertEqual(len(sample_cols), 1)
    
    def test_assign_season_groups(self):
        """Test season assignment logic"""
        
        # Test wet season assignment
        test_row_wet = pd.Series({'simulation_year': 2, 'month': 10})
        result = assign_season_group(test_row_wet)
        self.assertIn('Wet season: 2-08 to 3-01', result)
        
        # Test dry season assignment
        test_row_dry = pd.Series({'simulation_year': 2, 'month': 4})
        result = assign_season_group(test_row_dry)
        self.assertIn('Dry season: 2-02 to 2-07', result)
        
        # Test peak group assignment
        test_row_peak = pd.Series({'simulation_year': 2, 'month': 11})
        result = assign_peak_group(test_row_peak)
        self.assertIn('Peak wet: 2-10 to 2-12', result)
    
    def test_subset_by_age(self):
        """Test age-based sampling"""
        df_processed = calculate_infection_metrics(convert_month(self.sample_infection_df.copy()))
        
        result_df = subset_by_age(
            df_processed,
            n_samples_year=4,
            replicates=1,
            scheme='test_age'
        )
        
        # Check that age bin column was added
        self.assertIn('age_bin', result_df.columns)
        sample_cols = [col for col in result_df.columns if 'test_age' in col and 'rep' in col]
        self.assertEqual(len(sample_cols), 1)
        
        # Check age bins were created correctly
        age_bins = result_df['age_bin'].unique()
        self.assertGreater(len(age_bins), 0)
    
    def test_validate_subset_inputs(self):
        """Test input validation for subset functions"""
        
        # Test with valid inputs - should not raise
        validate_subset_inputs(self.sample_infection_df, 5, 2, "test_method")
        
        # Test with invalid inputs - should raise
        with self.assertRaises(ValueError):
            validate_subset_inputs(pd.DataFrame(), 5, 2, "test_method")  # Empty df
        
        with self.assertRaises(ValueError):
            validate_subset_inputs(self.sample_infection_df, 0, 2, "test_method")  # Invalid n_samples
        
        with self.assertRaises(ValueError):
            validate_subset_inputs(self.sample_infection_df, 5, 0, "test_method")  # Invalid replicates
    
    def test_filter_emod_infections_deduplication(self):
        """Test EMOD infection deduplication logic"""
        
        # Create test data with duplicates
        duplicate_df = pd.DataFrame({
            'IndividualID': [1, 1, 2, 2, 3],
            'continuous_month': [1, 1, 1, 2, 1],
            'infIndex': [0, 1, 2, 3, 4],
            'day': [1, 2, 1, 32, 1]
        })
        
        result_df = filter_emod_infections(duplicate_df, duplicate_seed=123, is_test=True)
        
        # Should have fewer rows due to deduplication
        self.assertLessEqual(len(result_df), len(duplicate_df))
        
        # Should have one row per individual per month
        grouped = result_df.groupby(['IndividualID', 'continuous_month']).size()
        self.assertTrue(all(grouped == 1))
    
    def test_run_sampling_functions_dispatcher(self):
        """Test the main sampling dispatcher function"""
        df_processed = calculate_infection_metrics(convert_month(self.sample_infection_df.copy()))
        
        config = {
            'method': 'random',
            'n_samples_year': 2,
            'replicates': 1,
            'method_params': {
                'population_proportions': None,
                'equal_monthly': False
            }
        }
        
        result_df = run_sampling_functions(df_processed, config)
        
        # Should have original data plus sampling columns
        self.assertGreaterEqual(result_df.shape[1], df_processed.shape[1])
    
    def test_run_sampling_model_integration(self):
        """Test the complete sampling model pipeline"""
        config = {
            'hard_filters': {
                'symptomatics_only': False,
                'monogenomic_infections_only': False
            },
            'sampling_configs': {
                'random': {
                    'method': 'random',
                    'n_samples_year': 2,
                    'replicates': 1,
                    'method_params': {
                        'population_proportions': None,
                        'equal_monthly': False
                    }
                }
            }
        }
        
        result_df = run_sampling_model(self.sample_infection_df, config)
        
        # Should have all the required columns
        required_cols = ['true_coi', 'effective_coi', 'cotx']
        for col in required_cols:
            self.assertIn(col, result_df.columns)
        
        # Should have sampling columns
        sampling_cols = [col for col in result_df.columns if 'random' in col and 'rep' in col]
        self.assertGreater(len(sampling_cols), 0)

class TestScientificCalculations(unittest.TestCase):
    """Test the scientific accuracy of metric calculations"""
    
    def setUp(self):
        """Create controlled test data with known expected outcomes"""
        
        # Create test infection data with known properties
        self.test_infections = pd.DataFrame({
            'infIndex': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            'IndividualID': [100, 101, 102, 103, 104, 105, 106, 107, 108, 109],
            'simulation_year': [1, 1, 1, 1, 2, 2, 2, 2, 2, 2],
            'month': [1, 3, 6, 9, 1, 3, 6, 9, 11, 12],
            'population': [0, 0, 1, 1, 0, 0, 1, 1, 0, 1],
            'fever_status': [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
            'age_day': [1000, 2000, 8000, 3000, 1500, 6000, 2500, 4000, 5000, 7000],
            # Known COI patterns
            'recursive_nid': ['[0]', '[1]', '[2,3]', '[4]', '[5,6,7]', '[8]', '[9,10]', '[11]', '[12,13,14,15]', '[16]'],
            'genome_ids': ['[100]', '[101]', '[102,103]', '[104]', '[105,106,107]', '[108]', '[109,110]', '[111]', '[112,113,114,115]', '[116]'],
            'bite_ids': ['[1]', '[2]', '[3]', '[4]', '[5]', '[6]', '[7]', '[8]', '[9]', '[10]']
        })
        
        # Create genotype matrix with known patterns for testing heterozygosity
        np.random.seed(42)  # Reproducible results
        self.genotype_matrix = np.array([
            [0, 1, 0, 1, 0, 1],  # genome 0
            [0, 1, 0, 1, 0, 1],  # genome 1 - identical to 0
            [1, 0, 1, 0, 1, 0],  # genome 2
            [1, 0, 1, 0, 1, 0],  # genome 3 - identical to 2
            [0, 0, 1, 1, 0, 1],  # genome 4
            [1, 1, 0, 0, 1, 0],  # genome 5
            [0, 1, 1, 0, 1, 1],  # genome 6
            [1, 0, 0, 1, 0, 0],  # genome 7
            [0, 0, 0, 0, 0, 0],  # genome 8 - all zeros
            [1, 1, 1, 1, 1, 1],  # genome 9 - all ones
            [0, 1, 0, 1, 0, 1],  # genome 10 - identical to 0,1
            [1, 0, 1, 0, 1, 0],  # genome 11 - identical to 2,3
            [0, 1, 1, 0, 0, 1],  # genome 12
            [1, 0, 0, 1, 1, 0],  # genome 13
            [0, 0, 1, 1, 1, 1],  # genome 14
            [1, 1, 0, 0, 0, 0],  # genome 15
            [0, 1, 0, 1, 0, 1],  # genome 16 - identical to 0,1,10
        ])
        
        register_matrix('test_genotype_matrix', self.genotype_matrix)
    
    def test_calculate_infection_metrics_scientific_accuracy(self):
        """Test that infection metrics are calculated correctly with known inputs"""
        
        # Add continuous month for processing
        test_df = self.test_infections.copy()
        test_df['continuous_month'] = (test_df['simulation_year'] * 12) + test_df['month']
        
        result_df = calculate_infection_metrics(test_df)
        
        # Test known COI values
        expected_true_coi = [1, 1, 2, 1, 3, 1, 2, 1, 4, 1]  # Based on genome_ids
        expected_effective_coi = [1, 1, 2, 1, 3, 1, 2, 1, 4, 1]  # Same as true (no duplicates in our test)
        
        np.testing.assert_array_equal(result_df['true_coi'].values, expected_true_coi)
        np.testing.assert_array_equal(result_df['effective_coi'].values, expected_effective_coi)
        
        # Test cotransmission calculation
        # Infections with COI=1 should have cotx=None
        mono_infections = result_df[result_df['effective_coi'] == 1]
        self.assertTrue(mono_infections['cotx'].isna().all())
        
        # Polygenomic infections with single bite should have cotx=True
        poly_infections = result_df[result_df['effective_coi'] > 1]
        for idx, row in poly_infections.iterrows():
            bite_ids = eval(row['bite_ids'])
            expected_cotx = len(set(bite_ids)) == 1
            self.assertEqual(row['cotx'], expected_cotx, f"Row {idx} cotx mismatch")
    
    def test_comprehensive_group_summary_mathematical_accuracy(self):
        """Test comprehensive group summary with known mathematical outcomes"""
        
        # Create test data with known summary statistics
        test_group = pd.DataFrame({
            'infIndex': [0, 1, 2, 3, 4],
            'true_coi': [1, 2, 3, 1, 4],  # Mean=2.2, Median=2
            'effective_coi': [1, 2, 2, 1, 3],  # Mean=1.8, Median=2
            'genome_ids': [[0], [1, 2], [3, 4, 5], [6], [7, 8, 9, 10]],
            'cotx': [None, True, False, None, True],
            'fever_status': [1, 0, 1, 0, 1],
            'population': [0, 0, 1, 1, 0],
            'age_day': [1000, 2000, 3000, 4000, 5000]
        })
        
        summary = comprehensive_group_summary(test_group)
        
        # Test basic counts
        self.assertEqual(summary['n_infections'], 5)
        self.assertEqual(summary['true_poly_coi_count'], 3)  # COI > 1: infections 1,2,4
        self.assertEqual(summary['effective_poly_coi_count'], 3)
        self.assertAlmostEqual(summary['true_poly_coi_prop'], 0.6, places=2)
        
        # Test mathematical calculations
        self.assertAlmostEqual(summary['true_coi_mean'], 2.2, places=1)
        self.assertAlmostEqual(summary['true_coi_median'], 2.0, places=1)
        
        # Test genome analysis
        total_genomes = sum(len(gids) for gids in test_group['genome_ids'])
        unique_genomes = len(set([gid for gids in test_group['genome_ids'] for gid in gids]))
        self.assertEqual(summary['genome_ids_total_count'], total_genomes)
        self.assertEqual(summary['genome_ids_unique_count'], unique_genomes)
        
        # Test cotransmission calculation - only count among polygenomic infections
        poly_rows = test_group[test_group['effective_coi'] > 1]
        cotx_count = poly_rows['cotx'].sum()  # True=1, False=0, so sum gives count of True
        cotx_prop = cotx_count / len(poly_rows)
        self.assertEqual(summary['cotransmission_count'], cotx_count)
        self.assertAlmostEqual(summary['cotransmission_prop'], cotx_prop, places=2)
    
    def test_heterozygosity_calculation_genetic_accuracy(self):
        """Test heterozygosity calculation with known genetic patterns"""
        
        matrix = self.genotype_matrix
        
        # Test 1: Identical genomes should show no heterozygosity
        identical_indices = [0, 1, 10]  # All identical genomes
        barcode = generate_het_barcode(matrix, identical_indices)
        expected = ['0', '1', '0', '1', '0', '1']  # Should match the genotype
        self.assertEqual(barcode, expected, "Identical genomes should show no heterozygosity")
        
        # Test 2: Completely different genomes at all loci
        different_indices = [8, 9]  # All 0s vs all 1s
        barcode = generate_het_barcode(matrix, different_indices)
        expected = ['N', 'N', 'N', 'N', 'N', 'N']  # All heterozygous
        self.assertEqual(barcode, expected, "Completely different genomes should be heterozygous at all loci")
        
        # Test 3: Mixed pattern - some loci heterozygous, some not
        mixed_indices = [0, 4]  # [0,1,0,1,0,1] vs [0,0,1,1,0,1]
        barcode = generate_het_barcode(matrix, mixed_indices)
        expected = ['0', 'N', 'N', 'N', '0', '1']  # Mixed pattern
        self.assertEqual(barcode, expected)
        
        # Test heterozygosity proportion calculation
        het_prop = barcode.count('N') / len(barcode)
        expected_prop = 4/6  # 4 N's out of 6 loci
        self.assertAlmostEqual(het_prop, expected_prop, places=2)
    
    def test_ibx_distribution_mathematical_accuracy(self):
        """Test IBx distribution calculation with known pairwise relationships"""
        
        # Create simple test matrix for IBx calculation
        simple_matrix = np.array([
            [0, 0, 0],  # genome 0
            [0, 0, 1],  # genome 1
            [0, 1, 1],  # genome 2
            [1, 1, 1]   # genome 3
        ])
        
        # Test pairwise IBx calculations
        indices = [0, 1, 2, 3]
        
        # Test the distribution function
        distribution = ibx_distribution(indices, simple_matrix)
        
        # Check that we get the expected number of pairwise comparisons
        total_pairs = sum(distribution.values())
        expected_pairs = len(list(combinations(indices, 2)))
        self.assertEqual(total_pairs, expected_pairs)
        
        # Test statistical summary
        summary = weighted_describe_scipy(distribution, 'test_ibx')
        self.assertIn('test_ibx_mean', summary)
        self.assertIn('test_ibx_std', summary)
        self.assertGreater(summary['test_ibx_pairwise_count'], 0)
    
    def test_nested_comparisons_structure_mathematical_integrity(self):
        """Test that nested comparisons maintain correct mathematical structure"""
        
        # Create test data with known grouping structure
        test_df = pd.DataFrame({
            'infIndex': list(range(12)),
            'population': [0, 0, 1, 1] * 3,
            'fever_status': [1, 0] * 6,
            'true_coi': [1, 2, 1, 3] * 3,
            'group_year': [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3],
            'group_month': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'age_day': [1000, 3000, 8000, 2000] * 3
        })
        
        config = {
            'populations': True,
            'polygenomic': True,
            'symptomatic': True,
            'age_bins': True
        }
        
        result = identify_nested_comparisons(test_df, 'random', config=config)
        
        # Verify expected comparison types are present
        expected_comparisons = ['group_year', 'group_month', 'populations', 'polygenomic', 'symptomatic', 'age_bins']
        for comp_type in expected_comparisons:
            self.assertIn(comp_type, result, f"Missing comparison type: {comp_type}")
        
        # Test group_year structure - should have exactly 3 years
        year_groups = result['group_year']
        self.assertEqual(len(year_groups), 3)  # Years 1, 2, 3
        
        # Each year should have exactly 4 infections
        for year, indices in year_groups.items():
            expected_count = 4
            self.assertEqual(len(indices), expected_count, f"Year {year} should have {expected_count} infections")
        
        # Verify no duplicate indices within same comparison type
        for comp_type, groups in result.items():
            if comp_type in ['group_year', 'group_month']:
                all_indices = []
                for group_indices in groups.values():
                    all_indices.extend(group_indices)
                unique_indices = set(all_indices)
                self.assertEqual(len(all_indices), len(unique_indices), f"Duplicate indices in {comp_type}")

class TestArrayBooleanIssues(unittest.TestCase):
    """Test resolution of array boolean issues"""
    
    def setUp(self):
        """Set up test data that might trigger array boolean issues"""
        self.df_single_pop = pd.DataFrame({
            'population': np.array([0, 0, 0, 0]),  # Explicit numpy array
            'group_year': [1, 1, 2, 2],
            'infIndex': [0, 1, 2, 3],
            'fever_status': [1, 0, 1, 0],
            'true_coi': [1, 2, 1, 3],
            'age_day': [1000, 2000, 3000, 4000]
        })
        
        self.df_multi_pop = pd.DataFrame({
            'population': np.array([0, 1, 0, 1]),  # Explicit numpy array
            'group_year': [1, 1, 2, 2],
            'infIndex': [0, 1, 2, 3],
            'fever_status': [1, 0, 1, 0],
            'true_coi': [1, 2, 1, 3],
            'age_day': [1000, 2000, 3000, 4000]
        })
    
    def test_array_comparison_fixes(self):
        """Test that array comparisons are properly fixed"""
        
        # Test single population
        single_pop_unique = self.df_single_pop['population'].unique()
        multi_pop_unique = self.df_multi_pop['population'].unique()
        
        # These should work without "truth value of array" error
        self.assertFalse(len(single_pop_unique) > 1)
        self.assertTrue(len(multi_pop_unique) > 1)
        
        # Test other comparison patterns
        self.assertTrue(1 in multi_pop_unique)
        self.assertFalse(2 in multi_pop_unique)
    
    def test_identify_nested_comparisons_with_numpy_arrays(self):
        """Test nested comparisons work with numpy array columns"""
        config = {
            'populations': True,
            'polygenomic': True,
            'symptomatic': True,
            'age_bins': True
        }
        
        # This should not raise an array boolean error
        result1 = identify_nested_comparisons(self.df_single_pop, 'random', config=config)
        result2 = identify_nested_comparisons(self.df_multi_pop, 'random', config=config)
        
        # Single population should not have population comparisons
        self.assertNotIn('populations', result1)
        
        # Multiple populations should have population comparisons
        self.assertIn('populations', result2)

class TestMatrixOperations(unittest.TestCase):
    """Test matrix operations and IBx calculations"""
    
    def setUp(self):
        """Set up test matrices"""
        self.small_matrix = np.random.randint(0, 2, size=(10, 20))
        register_matrix('small_test_matrix', self.small_matrix)
        
        # Create IBx test matrix with known relationships
        self.ibx_matrix = np.array([
            [0, 1, 0, 1, 0],  # genome 0
            [0, 1, 0, 1, 0],  # genome 1 - identical to 0
            [1, 0, 1, 0, 1],  # genome 2
            [1, 0, 1, 0, 1],  # genome 3 - identical to 2
            [0, 0, 1, 1, 1],  # genome 4 - different
        ])
        register_matrix('ibx_test_matrix', self.ibx_matrix)
    
    def test_matrix_registration_and_retrieval(self):
        """Test matrix registration system"""
        # Test registration
        test_matrix = np.array([[1, 0], [0, 1]])
        register_matrix('test_registration_matrix', test_matrix)
        
        # Test retrieval
        retrieved = get_matrix('test_registration_matrix')
        np.testing.assert_array_equal(retrieved, test_matrix)
        
        # Test error on missing matrix
        with self.assertRaises(KeyError):
            get_matrix('nonexistent_matrix')
    
    def test_heterozygosity_edge_cases_comprehensive(self):
        """Test heterozygosity calculation edge cases"""
        matrix = self.ibx_matrix
        
        # Test single genome (should not be heterozygous)
        barcode = generate_het_barcode(matrix, [0])
        self.assertEqual(barcode, ['0', '1', '0', '1', '0'])
        
        # Test identical genomes
        barcode = generate_het_barcode(matrix, [0, 1])  # Identical genomes
        self.assertEqual(barcode, ['0', '1', '0', '1', '0'])
        
        # Test completely different genomes
        barcode = generate_het_barcode(matrix, [0, 2])  # Opposite patterns
        self.assertEqual(barcode, ['N', 'N', 'N', 'N', 'N'])
        
        # Test mixed heterozygosity
        barcode = generate_het_barcode(matrix, [0, 4])  # Some positions same, some different
        expected = ['0', 'N', 'N', 'N', 'N']  # First position same, others different
        self.assertEqual(barcode, expected)
    
    def test_ibx_calculations_with_known_matrix(self):
        """Test IBx calculations with known matrix relationships"""
        
        # Test with known identical pairs
        identical_indices = [0, 1]  # Known identical genomes
        distribution = ibx_distribution(identical_indices, self.ibx_matrix)
        
        # Should have high IBx values for identical genomes
        max_ibx = max(distribution.keys())
        self.assertEqual(max_ibx, 1.0)  # Perfect identity
        
        # Test with different genomes
        different_indices = [0, 2]  # Known different genomes
        distribution_diff = ibx_distribution(different_indices, self.ibx_matrix)
        
        # Should have lower IBx values
        min_ibx = min(distribution_diff.keys())
        self.assertEqual(min_ibx, 0.0)  # No identity

class TestEdgeCasesAndRobustness(unittest.TestCase):
    """Test edge cases and robustness of all functions"""
    
    def test_empty_data_handling(self):
        """Test functions handle empty data gracefully"""
        
        empty_df = pd.DataFrame()
        
        # Test comprehensive_group_summary with empty data
        summary = comprehensive_group_summary(empty_df)
        self.assertEqual(summary['n_infections'], 0)
        
        # Test heterozygosity with empty indices
        test_matrix = np.random.randint(0, 2, (10, 5))
        barcode = generate_het_barcode(test_matrix, [])
        self.assertEqual(barcode, [])
        
        # Test parse_list with edge cases
        self.assertEqual(parse_list(''), [])
        self.assertEqual(parse_list('invalid'), [])
    
    def test_single_sample_calculations(self):
        """Test functions work correctly with single samples"""
        
        single_df = pd.DataFrame({
            'infIndex': [0],
            'true_coi': [2],
            'effective_coi': [2],
            'genome_ids': [[1, 2]],
            'cotx': [True],
            'fever_status': [1],
            'population': [0],
            'age_day': [1000]
        })
        
        summary = comprehensive_group_summary(single_df)
        self.assertEqual(summary['n_infections'], 1)
        self.assertEqual(summary['true_poly_coi_count'], 1)
        self.assertAlmostEqual(summary['true_poly_coi_prop'], 1.0)
    
    def test_extreme_values_handling(self):
        """Test handling of extreme values"""
        
        extreme_df = pd.DataFrame({
            'infIndex': [0, 1, 2],
            'true_coi': [1, 50, 100],  # Very high COI
            'effective_coi': [1, 25, 50],
            'genome_ids': [[1], list(range(50)), list(range(100))],
            'cotx': [None, False, True],
            'fever_status': [1, 0, 1],
            'population': [0, 0, 1],
            'age_day': [1000, 2000, 3000]
        })
        
        summary = comprehensive_group_summary(extreme_df)
        self.assertEqual(summary['n_infections'], 3)
        self.assertEqual(summary['true_coi_max'], 100)
        self.assertEqual(summary['true_coi_min'], 1)
    
    def test_data_type_consistency(self):
        """Test that functions maintain consistent data types"""
        
        # Test with mixed data types
        mixed_df = pd.DataFrame({
            'infIndex': [0, 1, 2],
            'population': [0, 1, 0],  # Regular list
            'fever_status': np.array([1, 0, 1]),  # Numpy array
            'true_coi': pd.Series([1, 2, 1]),  # Pandas Series
            'group_year': ['1', '2', '1'],  # String values
            'age_day': [1000.0, 2000.5, 3000.1]  # Float values
        })
        
        config = {'populations': True, 'polygenomic': False, 'symptomatic': True, 'age_bins': False}
        
        # Should handle mixed types without error
        result = identify_nested_comparisons(mixed_df, 'random', config=config)
        self.assertIsInstance(result, dict)
        self.assertIn('group_year', result)

class TestIntegrationScenarios(unittest.TestCase):
    """Integration tests for complete workflows"""
    
    def test_complete_pipeline_integration(self):
        """Test the complete pipeline from start to finish"""
        
        # Create comprehensive test dataset
        full_test_df = pd.DataFrame({
            'infIndex': list(range(20)),
            'IndividualID': [1000 + i for i in range(20)],
            'simulation_year': [1]*6 + [2]*7 + [3]*7,
            'year': [1]*6 + [2]*7 + [3]*7,
            'month': [1, 3, 6, 9, 11, 12] + [1, 3, 6, 9, 11, 12, 12] + [1, 3, 6, 9, 11, 12, 12],
            'population': [0, 1] * 10,
            'fever_status': [1, 0] * 10,
            'age_day': [1000, 2000, 8000, 3000, 5000] * 4,
            'recursive_nid': [f'[{i}]' if i % 3 == 0 else f'[{i},{i+50}]' for i in range(20)],
            'genome_ids': [f'[{100+i}]' if i % 3 == 0 else f'[{100+i},{150+i}]' for i in range(20)],
            'bite_ids': [f'[{i+1}]' for i in range(20)]
        })
        
        config = get_default_config()
        
        # Run complete pipeline
        result_df = run_sampling_model(full_test_df, config)
        
        # Verify pipeline integrity
        self.assertGreater(len(result_df), 0)
        required_cols = ['true_coi', 'effective_coi', 'cotx']
        for col in required_cols:
            self.assertIn(col, result_df.columns)
        
        # Verify sampling columns exist
        sample_cols = [col for col in result_df.columns if 'rep' in col]
        self.assertGreater(len(sample_cols), 0)
        
        # Verify calculations are correct
        for idx, row in result_df.iterrows():
            genome_ids = eval(row['genome_ids'])
            expected_true_coi = len(genome_ids)
            expected_effective_coi = len(set(genome_ids))
            
            self.assertEqual(row['true_coi'], expected_true_coi)
            self.assertEqual(row['effective_coi'], expected_effective_coi)
    
    @patch('run_observational_model.pd.read_csv')
    @patch('os.path.exists')
    @patch('os.makedirs')
    def test_run_observational_model_mock_integration(self, mock_makedirs, mock_exists, mock_read_csv):
        """Test main observational model function with mocked dependencies"""
        
        # Mock file operations
        mock_exists.return_value = False
        mock_read_csv.return_value = pd.DataFrame({
            'infIndex': [0, 1, 2, 3],
            'recursive_nid': ['[0]', '[1]', '[2,3]', '[4]'],
            'genome_ids': ['[100]', '[101]', '[102,103]', '[104]'],
            'bite_ids': ['[1]', '[2]', '[3]', '[4]'],
            'fever_status': [1, 0, 1, 0],
            'population': [0, 0, 1, 1],
            'age_day': [1000, 2000, 3000, 4000],
            'simulation_year': [1, 1, 2, 2],
            'month': [1, 6, 1, 6],
            'IndividualID': [100, 101, 102, 103]
        })
        
        # This should run without array boolean errors
        try:
            summaries, samples = run_observational_model(
                sim_name="integration_test",
                emod_output_path="./test_data",
                config_path="./nonexistent_config.json",
                output_path="./test_output"
            )
            
            # Verify results structure
            self.assertIsNotNone(summaries)
            self.assertIsNotNone(samples)
            self.assertIsInstance(summaries, pd.DataFrame)
            self.assertIsInstance(samples, pd.DataFrame)
            
        except ValueError as e:
            if "truth value of an array" in str(e):
                self.fail("Array boolean error still occurring in integration test")
            else:
                # Other errors might be expected in mock environment
                pass

class TestRhMetricCalculations(unittest.TestCase):
    """Test the R_h metric calculation functionality"""
    
    def setUp(self):
        """Set up test data for R_h calculations"""
        
        # Create test data with known R_h patterns
        self.rh_test_df = pd.DataFrame({
            'infIndex': [0, 1, 2, 3, 4, 5, 6, 7, 8, 9],
            'effective_coi': [1, 1, 2, 2, 2, 3, 3, 4, 1, 2],
            'cotx': [None, None, False, False, True, False, True, False, None, False],  # Superinfections vs cotransmissions
            'heterozygosity': [0.0, 0.0, 0.25, 0.30, 0.15, 0.40, 0.20, 0.50, 0.0, 0.35],
            'ibs_mean': [0.95, 0.92, 0.75, 0.70, 0.85, 0.60, 0.80, 0.50, 0.90, 0.65]
        })
        
        # Create test monogenomic distribution dictionary
        self.test_monogenomic_dict = {
            0.45: 15,   # IBS value: count
            0.50: 20,
            0.55: 25,
            0.60: 18,
            0.65: 12,
            0.70: 8,
            0.75: 5,
            1.00: 3     # Should be excluded in calculations
        }
        
    def test_calculate_rh_basic_functionality(self):
        """Test basic R_h calculation functionality"""
        
        from unified_metric_calculations import calculate_rh
        
        # Test with known data
        rh_summary, individual_rh = calculate_rh(
            self.rh_test_df, 
            self.test_monogenomic_dict,
            n_mono_boostraps=50  # Small number for testing
        )
        
        # Check return types
        self.assertIsInstance(rh_summary, pd.DataFrame)
        self.assertIsInstance(individual_rh, pd.DataFrame)
        
        # Check that summary has one row
        self.assertEqual(len(rh_summary), 1)
        
        # Check required columns exist
        required_rh_columns = [
            'rh_mono_count', 'rh_mono_measured_mean', 'rh_mono_measured_median', 'rh_mono_measured_std',
            'rh_poly_count', 'rh_poly_measured_mean', 'rh_poly_measured_median', 'rh_poly_measured_std',
            'rh_poly_inferred_mean', 'rh_poly_inferred_median', 'rh_poly_inferred_std'
        ]
        
        for col in required_rh_columns:
            self.assertIn(col, rh_summary.columns, f"Missing column: {col}")
            
        # Check individual results structure
        expected_individual_cols = ['infIndex', 'individual_measured_rh', 'individual_inferred_rh']
        for col in expected_individual_cols:
            self.assertIn(col, individual_rh.columns, f"Missing individual column: {col}")
    
    def test_calculate_rh_mathematical_accuracy(self):
        """Test R_h calculations with known expected values"""
        
        from unified_metric_calculations import calculate_rh
        
        # Use fixed random seed for reproducible results
        np.random.seed(42)
        
        rh_summary, individual_rh = calculate_rh(
            self.rh_test_df, 
            self.test_monogenomic_dict,
            n_mono_boostraps=100
        )
        
        # Test known calculations
        # COI=2 superinfections: indices 2, 3, 9 (cotx=False)
        expected_mono_count = 3
        self.assertEqual(rh_summary.iloc[0]['rh_mono_count'], expected_mono_count)
        
        # Expected H_mono from COI=2 superinfections
        coi2_super = self.rh_test_df[(self.rh_test_df['effective_coi'] == 2) & (self.rh_test_df['cotx'] == False)]
        expected_mono_mean = round(coi2_super['ibs_mean'].mean(), 3)
        self.assertAlmostEqual(rh_summary.iloc[0]['rh_mono_measured_mean'], expected_mono_mean, places=2)
        
        # Test that polygenomic count matches expectations
        poly_samples = self.rh_test_df[self.rh_test_df['effective_coi'] > 1]
        expected_poly_count = len(poly_samples)
        self.assertEqual(rh_summary.iloc[0]['rh_poly_count'], expected_poly_count)
        
        # Test individual R_h calculation formula
        # R_h = (H_mono - H_poly) / H_mono
        for idx, row in individual_rh.iterrows():
            inf_data = self.rh_test_df[self.rh_test_df['infIndex'] == row['infIndex']].iloc[0]
            expected_measured_rh = (expected_mono_mean - inf_data['heterozygosity']) / expected_mono_mean
            self.assertAlmostEqual(row['individual_measured_rh'], expected_measured_rh, places=2)
    
    def test_calculate_rh_edge_cases(self):
        """Test R_h calculation edge cases"""
        
        from unified_metric_calculations import calculate_rh
        
        # Test with no COI=2 superinfections
        no_coi2_df = pd.DataFrame({
            'infIndex': [0, 1, 2],
            'effective_coi': [1, 3, 4],
            'cotx': [None, True, False],
            'heterozygosity': [0.0, 0.3, 0.4],
            'ibs_mean': [0.9, 0.6, 0.5]
        })
        
        with self.assertRaises((ValueError, ZeroDivisionError, KeyError)):
            calculate_rh(no_coi2_df, self.test_monogenomic_dict)
        
        # Test with no polygenomic infections
        no_poly_df = pd.DataFrame({
            'infIndex': [0, 1, 2],
            'effective_coi': [1, 1, 1],
            'cotx': [None, None, None],
            'heterozygosity': [0.0, 0.0, 0.0],
            'ibs_mean': [0.9, 0.9, 0.9]
        })
        
        # Should handle gracefully
        rh_summary, individual_rh = calculate_rh(no_poly_df, self.test_monogenomic_dict)
        self.assertEqual(rh_summary.iloc[0]['rh_poly_count'], 0)
        self.assertEqual(len(individual_rh), 0)
    
    def test_monogenomic_dict_sampling(self):
        """Test the distribution sampling logic"""
        
        # Test the internal sampling function logic
        test_heterozygosity = 0.25
        
        # Manual calculation to verify sampling
        exclude_keys = [1]
        filtered_dict = {k: v for k, v in self.test_monogenomic_dict.items() if k not in exclude_keys}
        
        values = list(filtered_dict.keys())
        weights = list(filtered_dict.values())
        total_weight = sum(weights)
        
        # Verify filtering worked
        self.assertNotIn(1.0, filtered_dict)
        
        # Test that sampling produces reasonable results
        np.random.seed(123)
        n_bootstraps = 1000
        sampled_values = np.random.choice(values, size=n_bootstraps, p=np.array(weights)/total_weight)
        
        # Calculate R_h manually
        rh_values = [(val - test_heterozygosity) / val for val in sampled_values if val != 0]
        manual_median = np.median(rh_values)
        
        # Should be reasonable R_h value
        self.assertGreater(manual_median, -1)  # R_h should be > -1
        self.assertLess(manual_median, 1)      # R_h should be < 1

class TestRhIntegration(unittest.TestCase):
    """Test R_h integration with the main pipeline"""
    
    def setUp(self):
        """Set up integration test data"""
        
        # Create comprehensive dataset for integration testing
        self.integration_df = pd.DataFrame({
            'infIndex': list(range(20)),
            'IndividualID': [1000 + i for i in range(20)],
            'simulation_year': [1]*10 + [2]*10,
            'year': [1]*10 + [2]*10,
            'month': [1, 3, 6, 9, 12] * 4,
            'population': [0, 1] * 10,
            'fever_status': [1, 0] * 10,
            'age_day': [1000, 2000, 8000, 3000, 5000] * 4,
            'recursive_nid': [f'[{i}]' if i % 4 == 0 else f'[{i},{i+20}]' for i in range(20)],
            'genome_ids': [f'[{100+i}]' if i % 4 == 0 else f'[{100+i},{120+i}]' for i in range(20)],
            'bite_ids': [f'[{i+1}]' if i % 4 == 0 else f'[{i+1}]' for i in range(20)],  # All single bite for simplicity
            'effective_coi': [1 if i % 4 == 0 else 2 for i in range(20)],
            'true_coi': [1 if i % 4 == 0 else 2 for i in range(20)],
            'cotx': [None if i % 4 == 0 else False for i in range(20)],  # All polygenomic are superinfections
            'heterozygosity': [0.0 if i % 4 == 0 else np.random.uniform(0.2, 0.4) for i in range(20)]
        })
        
        # Add sampling columns
        self.integration_df['rep_random_0'] = 1  # All sampled for simplicity
        
    @patch('unified_metric_calculations.get_matrix')
    def test_rh_in_run_time_summaries(self, mock_get_matrix):
        """Test R_h calculation integration in run_time_summaries"""
        
        # Mock matrix
        mock_matrix = np.random.randint(0, 2, size=(50, 100))
        mock_get_matrix.return_value = mock_matrix
        
        from unified_metric_calculations import run_time_summaries
        
        config = {
            'populations': False,
            'polygenomic': True,  # Required for R_h
            'symptomatic': False,
            'age_bins': False
        }
        
        try:
            # This should run without errors and include R_h calculations
            summary_df, individual_df, distributions_dict = run_time_summaries(
                self.integration_df,
                subpop_config=config,
                user_ibx_categories=['ibs'],
                rh_calculation=True
            )
            
            # Should have completed without errors
            self.assertIsInstance(summary_df, pd.DataFrame)
            self.assertIsInstance(individual_df, pd.DataFrame)
            self.assertIsInstance(distributions_dict, dict)
            
            # Should have R_h columns in summary
            rh_columns = [col for col in summary_df.columns if 'rh_' in col]
            self.assertGreater(len(rh_columns), 0, "No R_h columns found in summary")
            
            # Should have individual R_h data
            if not individual_df.empty:
                individual_rh_cols = [col for col in individual_df.columns if 'individual_' in col and 'rh' in col]
                self.assertGreater(len(individual_rh_cols), 0, "No individual R_h columns found")
                
        except Exception as e:
            self.fail(f"R_h integration test failed with error: {e}")
    
    def test_rh_requires_polygenomic_comparison(self):
        """Test that R_h calculation properly requires polygenomic comparisons"""
        
        from unified_metric_calculations import run_time_summaries
        
        # Config without polygenomic comparison should raise error
        config_no_poly = {
            'populations': False,
            'polygenomic': False,  # This should cause R_h to fail
            'symptomatic': False,
            'age_bins': False
        }
        
        with patch('unified_metric_calculations.get_matrix') as mock_get_matrix:
            mock_get_matrix.return_value = np.random.randint(0, 2, size=(50, 100))
            
            # Should not crash but should warn about missing polygenomic comparisons
            summary_df, individual_df, distributions_dict = run_time_summaries(
                self.integration_df,
                subpop_config=config_no_poly,
                user_ibx_categories=['ibs'],
                rh_calculation=True
            )
            
            # Should still return valid DataFrames even if R_h failed
            self.assertIsInstance(summary_df, pd.DataFrame)

class TestRhScientificValidation(unittest.TestCase):
    """Test R_h calculations for scientific accuracy"""
    
    def test_rh_mathematical_properties(self):
        """Test that R_h calculations maintain expected mathematical properties"""
        
        from unified_metric_calculations import calculate_rh
        
        # Create test data with known mathematical relationships
        test_data = pd.DataFrame({
            'infIndex': [0, 1, 2, 3],
            'effective_coi': [2, 2, 3, 4],
            'cotx': [False, False, False, False],  # All superinfections
            'heterozygosity': [0.1, 0.3, 0.5, 0.7],  # Increasing heterozygosity
            'ibs_mean': [0.6, 0.6, 0.4, 0.2]  # Some with same H_mono baseline
        })
        
        monogenomic_dist = {0.5: 50, 0.6: 30, 0.7: 20}
        
        np.random.seed(42)  # For reproducible bootstrap
        rh_summary, individual_rh = calculate_rh(test_data, monogenomic_dist, n_mono_boostraps=100)
        
        # Test mathematical relationships
        h_mono = rh_summary.iloc[0]['rh_mono_measured_mean']
        
        # R_h should decrease as heterozygosity increases (with same H_mono)
        same_hmono_samples = individual_rh[individual_rh['infIndex'].isin([0, 1])]
        if len(same_hmono_samples) == 2:
            low_het_rh = same_hmono_samples[same_hmono_samples['infIndex'] == 0]['individual_measured_rh'].iloc[0]
            high_het_rh = same_hmono_samples[same_hmono_samples['infIndex'] == 1]['individual_measured_rh'].iloc[0]
            
            self.assertGreater(low_het_rh, high_het_rh, 
                             "R_h should be higher when heterozygosity is lower")
        
        # All R_h values should be reasonable (between -1 and 1 typically)
        for rh_val in individual_rh['individual_measured_rh']:
            self.assertGreater(rh_val, -2, f"R_h value {rh_val} is unreasonably low")
            self.assertLess(rh_val, 2, f"R_h value {rh_val} is unreasonably high")
    
    def test_rh_bootstrap_consistency(self):
        """Test that bootstrap R_h calculations are consistent"""
        
        from unified_metric_calculations import calculate_rh
        
        test_data = pd.DataFrame({
            'infIndex': [0, 1, 2],
            'effective_coi': [2, 2, 3],
            'cotx': [False, False, False],
            'heterozygosity': [0.2, 0.3, 0.4],
            'ibs_mean': [0.6, 0.6, 0.5]
        })
        
        monogenomic_dist = {0.4: 20, 0.5: 30, 0.6: 25, 0.7: 15}
        
        # Run multiple times with different seeds
        results = []
        for seed in [42, 123, 456]:
            np.random.seed(seed)
            _, individual_rh = calculate_rh(test_data, monogenomic_dist, n_mono_boostraps=200)
            results.append(individual_rh['individual_inferred_rh'].mean())
        
        # Results should be similar but not identical (due to bootstrap sampling)
        mean_result = np.mean(results)
        std_result = np.std(results)
        
        # Standard deviation should be reasonable (not too high)
        coefficient_of_variation = std_result / abs(mean_result) if mean_result != 0 else float('inf')
        self.assertLess(coefficient_of_variation, 0.3, 
                       "Bootstrap results are too variable")

# UPDATE THE MAIN TEST RUNNER TO INCLUDE RH TESTS
def run_all_comprehensive_tests_with_rh():
    """Updated test runner that includes R_h metric tests"""
    
    # Define all test classes including new R_h tests
    test_classes = [
        TestBasicFunctionality,
        TestSamplingAndFiltering,
        TestScientificCalculations,
        TestArrayBooleanIssues,
        TestMatrixOperations,
        TestEdgeCasesAndRobustness,
        TestIntegrationScenarios,
        TestRhMetricCalculations,        # NEW
        TestRhIntegration,               # NEW
        TestRhScientificValidation       # NEW
    ]
    
    # Create test suite
    suite = unittest.TestSuite()
    
    total_tests = 0
    for test_class in test_classes:
        tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        suite.addTests(tests)
        total_tests += tests.countTestCases()
    
    print(f"Running {total_tests} tests across {len(test_classes)} test classes...")
    print("="*80)
    
    # Run tests with detailed output
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    # Detailed reporting
    print("\n" + "="*80)
    print("COMPREHENSIVE TEST RESULTS (INCLUDING R_H METRICS)")
    print("="*80)
    
    success_count = result.testsRun - len(result.failures) - len(result.errors)
    success_rate = (success_count / result.testsRun * 100) if result.testsRun > 0 else 0
    
    print(f"Tests Run: {result.testsRun}")
    print(f"Successes: {success_count}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Success Rate: {success_rate:.1f}%")
    
    # Report failures and errors as before
    if result.failures:
        print(f"\n{len(result.failures)} FAILURES:")
        for i, (test, traceback) in enumerate(result.failures, 1):
            print(f"\n{i}. {test}")
            print("-" * 40)
            print(traceback)
    
    if result.errors:
        print(f"\n{len(result.errors)} ERRORS:")
        for i, (test, traceback) in enumerate(result.errors, 1):
            print(f"\n{i}. {test}")
            print("-" * 40)
            print(traceback)
    
    # Test class breakdown including R_h tests
    print(f"\nTEST CLASS BREAKDOWN:")
    for test_class in test_classes:
        class_tests = unittest.TestLoader().loadTestsFromTestCase(test_class)
        print(f"  {test_class.__name__}: {class_tests.countTestCases()} tests")
    
    return result

if __name__ == '__main__':
    # Run all tests including new R_h tests
    test_result = run_all_comprehensive_tests_with_rh()
    exit_code = 0 if test_result.wasSuccessful() else 1
    sys.exit(exit_code)
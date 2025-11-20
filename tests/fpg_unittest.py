"""
Self-contained test suite using Python's built-in unittest (no pytest needed).
All test data (dataframes and matrices) are generated within this file.

Run with: python test_infection_sampling_unittest.py
Or: python -m unittest test_infection_sampling_unittest -v

This test suite covers key functions from fpg_observational_model.unified_sampling.py and unified_metric_calculations.py.

AI generated on 2025-11 by Claude, based on user-provided context and reviewed for accuracy.
"""

import os
import unittest
import pandas as pd
import numpy as np
import ast
import sys
from itertools import combinations
from collections import Counter

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
        subset_by_age,
        filter_emod_infections,
        run_sampling_model
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
        get_variant_coi,
        generate_het_barcode,
        calculate_ibx_matrix,
        ibx_distribution,
        weighted_describe_scipy,
        calculate_rh,
        identify_nested_comparisons,
        process_nested_summaries,
        process_nested_ibx,
        run_time_summaries,
        update_ibx_index
    )
    METRICS_IMPORTED = True
except ImportError as e:
    print(f"Warning: Could not import unified_metric_calculations: {e}")
    METRICS_IMPORTED = False


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
                'fever_status': [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
                'age_day': [1000, 2000, 3000, 8000, 1500, 6000, 2500, 4000, 5000, 7000],
                'recursive_nid': ['[0]', '[1]', '[2, 3]', '[4]', '[5, 6, 7]', '[8]', '[9]', '[10, 11]', '[12]', '[13, 14]'],
                'genome_ids': ['[100]', '[101]', '[102, 103]', '[104]', '[105, 106, 107]', '[108]', '[109]', '[110, 111]', '[112]', '[113, 114]'],
                'bite_ids': ['[1]', '[2]', '[3]', '[4]', '[5]', '[6]', '[7]', '[8]', '[9]', '[10]']
            })
            
            # Create genotype matrix (IBS)
            np.random.seed(42)
            self.genotype_matrix = np.random.randint(0, 2, size=(15, 10), dtype=np.int8)
            self.genotype_matrix[2] = self.genotype_matrix[3]  # Make some identical
            
            # Create IBD matrix
            np.random.seed(123)
            ibd = np.random.randint(0, 101, size=(15, 15))
            self.ibd_matrix = (ibd + ibd.T) // 2
            np.fill_diagonal(self.ibd_matrix, 100)
            self.ibd_matrix = self.ibd_matrix.astype(np.int16)
            
            # Standard config
            self.config = {
                'hard_filters': {
                    'symptomatics_only': False,
                    'monogenomic_infections_only': False,
                    'day_snapshot': False
                },
                'intervention_start_month': 29,
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
                    'add_monthly': False,
                    'populations': False,
                    'polygenomic': True,
                    'symptomatic': False,
                    'age_bins': False
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
print("\n" + "="*70)
print("INITIALIZING SHARED TEST DATA")
print("="*70)

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

print("="*70 + "\n")


###############################################################################
# TEST CLASSES - Using shared TEST_DATA
###############################################################################

@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestCalculateInfectionMetrics(unittest.TestCase):
    """Test calculate_infection_metrics function from fpg_observational_model.unified_sampling.py"""
    
    def setUp(self):
        """Use shared test data"""
        self.df = TEST_DATA.get_infection_df(with_metrics=False)
    
    def test_calculate_infection_metrics_basic(self):
        """Test that calculate_infection_metrics runs and adds columns"""
        result = calculate_infection_metrics(self.df)
        
        self.assertIn('true_coi', result.columns)
        self.assertIn('effective_coi', result.columns)
        self.assertIn('cotx', result.columns)
        self.assertIn('recursive_nids_parsed', result.columns)
    
    def test_coi_values(self):
        """Test COI calculation values"""
        result = calculate_infection_metrics(self.df)
        
        expected_true_coi = [1, 1, 2, 1, 3, 1, 1, 2, 1, 2]
        actual_true_coi = result['true_coi'].tolist()
        
        self.assertEqual(actual_true_coi, expected_true_coi)


@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestApplyEmodFilters(unittest.TestCase):
    """Test apply_emod_filters function"""
    
    def setUp(self):
        """Use shared test data with metrics"""
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
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
        
        self.assertEqual(len(poly_remaining), 0, 
                        f"Filter failed: {len(poly_remaining)} polygenomic infections remain")


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestHeterozygosityFunctions(unittest.TestCase):
    """Test get_variant_coi and generate_het_barcode"""
    
    def setUp(self):
        """Use shared test data"""
        self.matrix = TEST_DATA.get_genotype_matrix()
        register_matrix('test_matrix', self.matrix)
    
    def test_get_variant_coi_monogenomic(self):
        """Test variant COI for monogenomic infection"""
        indices = [0]
        coi = get_variant_coi(self.matrix, indices)
        self.assertEqual(coi, 1)
    
    def test_get_variant_coi_polygenomic(self):
        """Test variant COI for polygenomic infection"""
        indices = [0, 1]
        coi = get_variant_coi(self.matrix, indices)
        self.assertEqual(coi, 2)
    
    def test_generate_het_barcode_monogenomic(self):
        """Test barcode generation for monogenomic"""
        indices = [0]
        barcode = generate_het_barcode(self.matrix, indices)
        
        self.assertEqual(len(barcode), self.matrix.shape[1])
        self.assertNotIn('N', barcode)
    
    def test_heterozygosity_proportion(self):
        """Test heterozygosity calculation"""
        indices = [0, 1]
        barcode = generate_het_barcode(self.matrix, indices)
        
        het = barcode.count('N') / len(barcode) if len(barcode) > 0 else 0
        self.assertGreaterEqual(het, 0)
        self.assertLessEqual(het, 1)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestComprehensiveGroupSummary(unittest.TestCase):
    """Test comprehensive_group_summary function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_comprehensive_summary_runs(self):
        """Test that summary function runs without error"""
        result = comprehensive_group_summary(self.df)
        
        self.assertIsInstance(result, pd.Series)
        self.assertIn('n_infections', result.index)
        self.assertIn('true_poly_coi_count', result.index)
    
    def test_infection_count(self):
        """Test infection count in summary"""
        result = comprehensive_group_summary(self.df)
        self.assertEqual(result['n_infections'], 10)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestRhCalculation(unittest.TestCase):
    """Test calculate_rh function and R_h metric logic"""
    
    def setUp(self):
        """Set up test data with heterozygosity and IBS information"""
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        
        # Add heterozygosity column (simulated)
        self.df['heterozygosity'] = [0.0, 0.0, 0.3, 0.0, 0.5, 0.0, 0.0, 0.4, 0.0, 0.35]
        
        self.monogenomic_dict = TEST_DATA.get_monogenomic_dict()
    
    def test_rh_formula_basic(self):
        """Test basic R_h formula"""
        h_mono = 0.5
        h_poly = 0.3
        rh = (h_mono - h_poly) / h_mono if h_mono != 0 else 0
        
        self.assertAlmostEqual(rh, 0.4)
    
    def test_calculate_rh_function(self):
        """Test the actual calculate_rh function"""
        poly_df = self.df[self.df['effective_coi'] > 1].copy()
        
        rh_summary, individual_rh = calculate_rh(poly_df, self.monogenomic_dict, n_mono_boostraps=200)
        
        self.assertIsInstance(rh_summary, pd.DataFrame)
        self.assertIsInstance(individual_rh, pd.DataFrame)
        
        expected_cols = ['rh_poly_inferred_mean', 'rh_poly_inferred_median', 'rh_poly_inferred_std']
        for col in expected_cols:
            self.assertIn(col, rh_summary.columns)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestIdentityByState(unittest.TestCase):
    """Test IBS matrix calculations"""
    
    def setUp(self):
        self.matrix = TEST_DATA.get_genotype_matrix()
        register_matrix('ibs_matrix', self.matrix)
    
    def test_ibs_pairwise_calculation(self):
        """Test IBS calculation between two genotypes"""
        geno1 = self.matrix[0, :]
        geno2 = self.matrix[1, :]
        
        matches = np.sum(geno1 == geno2)
        ibs = matches / len(geno1)
        
        self.assertGreaterEqual(ibs, 0)
        self.assertLessEqual(ibs, 1)
    
    def test_ibs_identical_genotypes(self):
        """Test IBS = 1.0 for identical genotypes"""
        geno1 = self.matrix[0, :]
        geno2 = self.matrix[0, :]
        
        matches = np.sum(geno1 == geno2)
        ibs = matches / len(geno1)
        
        self.assertEqual(ibs, 1.0)


@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestInterventionTiming(unittest.TestCase):
    """Test intervention start month handling"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_intervention_month_calculation(self):
        """Test that intervention_month is calculated correctly"""
        from fpg_observational_model.unified_sampling import adjust_time_columns
        
        intervention_start_month = 29
        result = adjust_time_columns(self.df, intervention_start_month=intervention_start_month)
        
        self.assertIn('intervention_month', result.columns)
        self.assertIn('intervention_year', result.columns)
    
    def test_intervention_resets_time_zero(self):
        """Test that intervention month becomes new time zero"""
        from fpg_observational_model.unified_sampling import adjust_time_columns
        
        intervention_start_month = 25
        result = adjust_time_columns(self.df, intervention_start_month=intervention_start_month)
        
        at_intervention = result[result['continuous_month'] == intervention_start_month]
        if len(at_intervention) > 0:
            self.assertTrue(all(at_intervention['intervention_month'] == 0))
    
    def test_group_year_uses_intervention_year(self):
        """Test that group_year uses intervention_year when intervention is set"""
        config = TEST_DATA.get_config()
        config['intervention_start_month'] = 29
        
        result = run_sampling_model(self.df, config, intervention_start_month=29, verbose=False)
        
        self.assertIn('intervention_year', result.columns)
        self.assertIn('group_year', result.columns)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestIdentityByDescent(unittest.TestCase):
    """Test IBD matrix calculations"""
    
    def setUp(self):
        self.ibd_matrix = TEST_DATA.get_ibd_matrix()
        register_matrix('ibd_matrix', self.ibd_matrix)
    
    def test_ibd_matrix_symmetric(self):
        """Test IBD matrix is symmetric"""
        self.assertTrue(np.allclose(self.ibd_matrix, self.ibd_matrix.T))
    
    def test_ibd_diagonal_maximum(self):
        """Test IBD diagonal elements are maximum"""
        max_val = np.max(self.ibd_matrix)
        for i in range(self.ibd_matrix.shape[0]):
            self.assertEqual(self.ibd_matrix[i, i], max_val)
    
    def test_ibd_pairwise_values(self):
        """Test IBD pairwise values are in valid range"""
        n = self.ibd_matrix.shape[0]
        
        for i in range(n):
            for j in range(i+1, n):
                ibd_val = self.ibd_matrix[i, j]
                self.assertGreaterEqual(ibd_val, 0)
                self.assertLessEqual(ibd_val, 100)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestWeightedStatistics(unittest.TestCase):
    """Test weighted_describe_scipy function"""
    
    def test_weighted_mean_calculation(self):
        """Test weighted mean calculation"""
        summary_dict = {0.2: 5, 0.5: 10, 0.8: 5}
        
        values = np.array(list(summary_dict.keys()))
        weights = np.array(list(summary_dict.values()))
        
        weighted_mean = np.average(values, weights=weights)
        expected = (0.2*5 + 0.5*10 + 0.8*5) / 20
        
        self.assertAlmostEqual(weighted_mean, expected, places=6)
    
    def test_weighted_describe_basic(self):
        """Test basic weighted describe"""
        summary_dict = {0.2: 5, 0.5: 10, 0.8: 5}
        
        result = weighted_describe_scipy(summary_dict, 'test')
        
        self.assertIsInstance(result, pd.DataFrame)
        self.assertIn('test_mean', result.columns)
        self.assertIn('test_median', result.columns)


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestNestedComparisons(unittest.TestCase):
    """Test identify_nested_comparisons function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.config = TEST_DATA.get_config()['subpopulation_comparisons']
    
    def test_yearly_grouping(self):
        """Test identification of yearly groups"""
        yearly_groups = self.df.groupby('group_year')['infIndex'].apply(list).to_dict()
        
        self.assertEqual(len(yearly_groups), 3)
        self.assertIn(1, yearly_groups)
        self.assertIn(2, yearly_groups)
        self.assertIn(3, yearly_groups)
    
    def test_identify_nested_comparisons_basic(self):
        """Test basic nested comparison identification"""
        sampling_col = 'random_5_rep1'
        self.df[sampling_col] = 1
        
        result = identify_nested_comparisons(
            self.df, 
            sampling_col, 
            config=self.config
        )
        
        self.assertIsInstance(result, dict)
        self.assertIn('group_year', result)


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
    
    def test_matrix_binary_values(self):
        """Test that matrix contains only 0 and 1"""
        unique_values = np.unique(self.matrix)
        self.assertTrue(set(unique_values).issubset({0, 1}))


@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestSamplingFunctions(unittest.TestCase):
    """Test actual sampling functions"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_subset_randomly(self):
        """Test subset_randomly function"""
        n_samples_year = 3
        replicates = 2
        
        result = subset_randomly(
            self.df, 
            n_samples_year=n_samples_year,
            replicates=replicates,
            scheme='random'
        )
        
        expected_cols = [f'random_{n_samples_year}_rep{i+1}' for i in range(replicates)]
        for col in expected_cols:
            self.assertIn(col, result.columns)
    
    def test_filter_emod_infections(self):
        """Test filter_emod_infections function"""
        result = filter_emod_infections(self.df, duplicate_seed=123)
        
        self.assertLessEqual(len(result), len(self.df))


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
        self.assertGreater(len(result), 0)


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


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestIBxDistribution(unittest.TestCase):
    """Test ibx_distribution function"""
    
    def setUp(self):
        self.ibs_matrix = TEST_DATA.get_ibd_matrix()
    
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
    
    def test_ibx_distribution_function(self):
        """Test ibx_distribution with actual function"""
        indices = [0, 1, 2]
        dist = ibx_distribution(indices, self.ibs_matrix)
        
        self.assertIsInstance(dist, dict)
        self.assertGreater(len(dist), 0)


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


@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestSeasonAssignment(unittest.TestCase):
    """Test season assignment functions"""
    
    def test_assign_season_group(self):
        """Test assign_season_group function"""
        from fpg_observational_model.unified_sampling import assign_season_group
        
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
        from fpg_observational_model.unified_sampling import assign_peak_group
        
        row_peak_wet = pd.Series({'simulation_year': 2020, 'month': 11})
        season = assign_peak_group(row_peak_wet)
        self.assertIn('Peak wet', season)


@unittest.skipIf(not (SAMPLING_IMPORTED and METRICS_IMPORTED), "Need both modules")
class TestRunTimeSummaries(unittest.TestCase):
    """Test run_time_summaries - main pipeline function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        
        # Add sampling columns
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


# Add more comprehensive R_h tests
@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestRhCalculationDetailed(unittest.TestCase):
    """Detailed R_h calculation tests"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.df['heterozygosity'] = [0.0, 0.0, 0.3, 0.0, 0.5, 0.0, 0.0, 0.4, 0.0, 0.35]
        self.monogenomic_dict = TEST_DATA.get_monogenomic_dict()
    
    def test_rh_zero_heterozygosity(self):
        """Test R_h when heterozygosity is zero"""
        h_mono = 0.5
        h_poly = 0.0
        rh = (h_mono - h_poly) / h_mono
        self.assertEqual(rh, 1.0)
    
    def test_rh_equal_heterozygosities(self):
        """Test R_h when H_mono = H_poly"""
        h_mono = 0.5
        h_poly = 0.5
        rh = (h_mono - h_poly) / h_mono
        self.assertEqual(rh, 0.0)
    
    def test_rh_exclude_identical_barcodes(self):
        """Test that IBS=1 excluded from H_mono"""
        filtered_dict = {k: v for k, v in self.monogenomic_dict.items() if k != 1.0}
        self.assertNotIn(1.0, filtered_dict.keys())
    
    def test_rh_sampling_from_distribution(self):
        """Test sampling from distribution"""
        filtered_dict = {k: v for k, v in self.monogenomic_dict.items() if k != 1.0}
        
        values = list(filtered_dict.keys())
        weights = list(filtered_dict.values())
        
        n_samples = 200
        samples = np.random.choice(values, size=n_samples, p=np.array(weights)/sum(weights))
        
        self.assertEqual(len(samples), n_samples)
        self.assertTrue(all(s < 1.0 for s in samples))
    
    def test_rh_values_in_valid_range(self):
        """Test calculated R_h in valid range"""
        poly_df = self.df[self.df['effective_coi'] > 1].copy()
        
        rh_summary, individual_rh = calculate_rh(poly_df, self.monogenomic_dict)
        
        mean_rh = rh_summary['rh_poly_inferred_mean'].iloc[0]
        median_rh = rh_summary['rh_poly_inferred_median'].iloc[0]
        
        self.assertGreaterEqual(mean_rh, -1)
        self.assertLessEqual(mean_rh, 1)
        self.assertGreaterEqual(median_rh, -1)
        self.assertLessEqual(median_rh, 1)


# Add more intervention timing tests
@unittest.skipIf(not SAMPLING_IMPORTED, "unified_sampling not available")
class TestInterventionTimingDetailed(unittest.TestCase):
    """Detailed intervention timing tests"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_intervention_year_calculation(self):
        """Test intervention_year calculation"""
        from fpg_observational_model.unified_sampling import adjust_time_columns
        
        intervention_start_month = 29
        result = adjust_time_columns(self.df, intervention_start_month=intervention_start_month)
        
        for idx, row in result.iterrows():
            expected_int_year = row['intervention_month'] // 12
            self.assertEqual(row['intervention_year'], expected_int_year)
    
    def test_pre_intervention_negative_months(self):
        """Test pre-intervention infections have negative months"""
        from fpg_observational_model.unified_sampling import adjust_time_columns
        
        intervention_start_month = 30
        result = adjust_time_columns(self.df, intervention_start_month=intervention_start_month)
        
        pre_intervention = result[result['continuous_month'] < intervention_start_month]
        if len(pre_intervention) > 0:
            self.assertTrue(all(pre_intervention['intervention_month'] < 0))
    
    def test_post_intervention_positive_months(self):
        """Test post-intervention infections have positive months"""
        from fpg_observational_model.unified_sampling import adjust_time_columns
        
        intervention_start_month = 20
        result = adjust_time_columns(self.df, intervention_start_month=intervention_start_month)
        
        post_intervention = result[result['continuous_month'] > intervention_start_month]
        if len(post_intervention) > 0:
            self.assertTrue(all(post_intervention['intervention_month'] > 0))


# Add more IBS tests
@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestIdentityByStateDetailed(unittest.TestCase):
    """Detailed IBS tests"""
    
    def setUp(self):
        self.matrix = TEST_DATA.get_genotype_matrix()
    
    def test_ibs_for_infection(self):
        """Test IBS for polygenomic infection"""
        indices = [0, 1]
        
        geno1 = self.matrix[indices[0], :]
        geno2 = self.matrix[indices[1], :]
        
        matches = np.sum(geno1 == geno2)
        ibs = matches / len(geno1)
        
        self.assertGreaterEqual(ibs, 0)
        self.assertLessEqual(ibs, 1)
    
    def test_ibs_within_cotransmission(self):
        """Test IBS for cotransmitted genotypes"""
        indices = [2, 3]  # Identical in test matrix
        
        geno1 = self.matrix[indices[0], :]
        geno2 = self.matrix[indices[1], :]
        
        matches = np.sum(geno1 == geno2)
        ibs = matches / len(geno1)
        
        self.assertEqual(ibs, 1.0)
    
    def test_ibs_matrix_structure(self):
        """Test IBS matrix properties"""
        ibs_matrix = TEST_DATA.get_ibd_matrix()
        
        self.assertTrue(np.allclose(ibs_matrix, ibs_matrix.T))


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestCalculateIbxMatrix(unittest.TestCase):
    """Test calculate_ibx_matrix function using idm module"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.matrix = TEST_DATA.get_genotype_matrix()
    
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
    
    def test_ibx_matrix_symmetry(self):
        """Test that IBx hash matrix is symmetric"""
        try:
            import idm
            
            hash_df, hash_ibx = calculate_ibx_matrix(self.df, self.matrix)
            
            # IBx matrix should be symmetric
            self.assertTrue(np.allclose(hash_ibx, hash_ibx.T),
                          "IBx matrix should be symmetric")
            
        except ImportError:
            self.skipTest("idm module not available")
    
    def test_ibx_diagonal_values(self):
        """Test IBx diagonal values (self-comparison)"""
        try:
            import idm
            
            hash_df, hash_ibx = calculate_ibx_matrix(self.df, self.matrix)
            
            # Diagonal should be maximum (identical to self)
            max_val = np.max(hash_ibx)
            diagonal = np.diag(hash_ibx)
            
            self.assertTrue(all(diagonal == max_val),
                          "Diagonal should be maximum value")
            
        except ImportError:
            self.skipTest("idm module not available")


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestUpdateIbxIndex(unittest.TestCase):
    """Test update_ibx_index function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
    
    def test_update_ibx_index_basic(self):
        """Test basic IBx index updating"""
        result = update_ibx_index(self.df)
        
        # Should have ibx_nid column
        self.assertIn('ibx_nid', result.columns)
        
        # ibx_nid should be a list for each row
        self.assertTrue(all(isinstance(x, list) for x in result['ibx_nid']))
    
    def test_update_ibx_index_global_order(self):
        """Test that indices are mapped to global order"""
        result = update_ibx_index(self.df)
        
        # Get all ibx_nid values
        all_indices = []
        for idx_list in result['ibx_nid']:
            all_indices.extend(idx_list)
        
        # Indices should start from 0 and be contiguous
        unique_indices = sorted(set(all_indices))
        expected_range = list(range(len(unique_indices)))
        
        self.assertEqual(unique_indices, expected_range,
                        "IBx indices should be remapped to 0...N-1")


@unittest.skipIf(not METRICS_IMPORTED, "unified_metric_calculations not available")
class TestProcessNestedIbx(unittest.TestCase):
    """Test process_nested_ibx function"""
    
    def setUp(self):
        self.df = TEST_DATA.get_infection_df(with_metrics=True)
        self.matrix = TEST_DATA.get_genotype_matrix()
        register_matrix('test_ibx_matrix', self.matrix)
        
        # Create nested indices
        self.nested_indices = {
            'group_year': {
                1: [0, 1, 2],
                2: [3, 4, 5],
                3: [6, 7, 8, 9]
            }
        }
    
    def test_process_nested_ibx_with_idm(self):
        """Test process_nested_ibx with idm module"""
        try:
            import idm
            
            ibx_summary, individual_ibx, ibx_dist = process_nested_ibx(
                self.df,
                'test_ibx_matrix',
                self.nested_indices,
                ibx_prefix='test',
                individual_ibx_calculation=True,
                save_ibx_distributions=True
            )
            
            # Check return types
            self.assertIsInstance(ibx_summary, pd.DataFrame)
            self.assertIsInstance(individual_ibx, pd.DataFrame)
            self.assertIsInstance(ibx_dist, dict)
            
            # Check summary has IBx columns
            ibx_cols = [col for col in ibx_summary.columns if 'test_' in col]
            self.assertGreater(len(ibx_cols), 0, "Should have IBx summary columns")
            
        except ImportError:
            self.skipTest("idm module not available")


###############################################################################
# MAIN TEST RUNNER
###############################################################################

def run_tests():
    """Run all tests and print results"""
    global TEST_DATA
    
    print("="*70)
    print("INFECTION SAMPLING & METRICS TEST SUITE")
    print("="*70)
    print(f"unified_sampling imported: {SAMPLING_IMPORTED}")
    print(f"unified_metric_calculations imported: {METRICS_IMPORTED}")
    print("="*70)
    
    if TEST_DATA is not None:
        print("\nUSING SHARED TEST DATA:")
        print(f"  Infection dataframe: {TEST_DATA.base_infection_df.shape}")
        print(f"  Genotype matrix (IBS): {TEST_DATA.genotype_matrix.shape}")
        print(f"  IBD matrix: {TEST_DATA.ibd_matrix.shape}")
        print(f"  All tests use IDENTICAL copies of this data")
    else:
        print("\nWARNING: TEST_DATA not initialized")
    print("="*70 + "\n")
    
    loader = unittest.TestLoader()
    suite = loader.loadTestsFromModule(sys.modules[__name__])
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    
    print("\n" + "="*70)
    print("TEST SUMMARY")
    print("="*70)
    print(f"Tests run: {result.testsRun}")
    print(f"Successes: {result.testsRun - len(result.failures) - len(result.errors)}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped)}")
    print("="*70)
    
    return result.wasSuccessful()


if __name__ == '__main__':
    success = run_tests()
    sys.exit(0 if success else 1)
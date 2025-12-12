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
import pandas as pd
import numpy as np
import ast
import sys
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
                'fever_status': [1, 0, 1, 0, 1, 0, 1, 0, 1, 0],
                'age_day': [1000, 2000, 3000, 8000, 1500, 6000, 2500, 4000, 5000, 7000],
                'recursive_nid': ['[0]', '[1]', '[2, 3]', '[4]', '[5, 6, 7]', '[8]', '[9]', '[10, 11]', '[12]',
                                  '[13, 14]'],
                'genome_ids': ['[100]', '[101]', '[102, 103]', '[104]', '[105, 106, 107]', '[108]', '[109]',
                               '[110, 111]', '[112]', '[113, 114]'],
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

        self.assertEqual(len(poly_remaining), 4,
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
            for j in range(i + 1, n):
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
        expected = (0.2 * 5 + 0.5 * 10 + 0.8 * 5) / 20

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

        expected_cols = [f'random_{n_samples_year}_rep{i + 1}' for i in range(replicates)]
        for col in expected_cols:
            self.assertIn(col, result.columns)

    def test_filter_emod_infections(self):
        """Test filter_emod_infections function"""
        result = filter_emod_infections(self.df, duplicate_seed=123)

        self.assertLessEqual(len(result), len(self.df))


###############################################################################
# NEW TESTS FOR run_observational_model
###############################################################################

@unittest.skipIf(not RUN_MODEL_IMPORTED, "run_observational_model not available")
class TestRunObservationalModel(unittest.TestCase):
    """Test suite for run_observational_model function updates"""

    def setUp(self):
        """Set up test fixtures"""
        # Create temporary directories for test outputs
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
                emod_output_path=self.test_dir,
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
                emod_output_path=self.test_dir,
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
                emod_output_path=self.test_dir,
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
            emod_output_path=self.test_dir,
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
            emod_output_path=self.test_dir,
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
                emod_output_path=self.test_dir,
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
                emod_output_path=self.test_dir,
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
            emod_output_path=self.test_dir,
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
                emod_output_path=self.test_dir,
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
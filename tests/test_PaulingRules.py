from PaulingRules import Pauling1, Pauling2, PaulingConnection, Pauling3and4, Pauling3, Pauling4, Pauling5, \
    RuleCannotBeAnalyzedError, is_an_oxide_and_no_env_for_O, FrequencyEnvironmentPauling1, Pauling0, \
    get_entropy_from_frequencies, get_most_frequent_environment, get_mean_CN_from_frequencies, \
    Pauling2_optimized_environments
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.bond_valence import BVAnalyzer
from collections import Counter
import unittest
import json
import math
import tempfile
import os



class is_an_oxide_and_no_env_for_O_Test(unittest.TestCase):

    def setUp(self):
        with open("mp-7000.json", "r") as f:
            dict_lse = json.load(f)
        lse = LightStructureEnvironments.from_dict(dict_lse)
        struct = lse.structure
        bva = BVAnalyzer()
        valences = bva.get_valences(structure=struct)
        lgf = LocalGeometryFinder()
        lgf.setup_structure(structure=struct)
        se = lgf.compute_structure_environments(maximum_distance_factor=1.41, only_cations=False, valences=valences)
        strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
        self.lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)

        with open("mp-5634.json", "r") as f:
            dict_lse2 = json.load(f)
        self.lse2 = LightStructureEnvironments.from_dict(dict_lse2)

    def test_is_an_oxide_and_no_env_for_O(self):
        with self.assertRaises(ValueError):
            is_an_oxide_and_no_env_for_O(self.lse)
        with self.assertRaises(ValueError):
            is_an_oxide_and_no_env_for_O(self.lse2)


class Get_entropy_from_frequencies_Test(unittest.TestCase):
    def test_get_entropy_from_frequencies(self):
        compare1 = 1 - (-math.log(0.5, 2) * 0.5 * 2 / (-math.log(1 / 66.0, 2) * 1.0 / 66.0 * 66.0))
        self.assertAlmostEqual(get_entropy_from_frequencies({'Ga': ['O:6', 'O:6', 'T:4', 'T:4']})["Ga"], compare1)
        self.assertDictEqual(get_entropy_from_frequencies({'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}),
                             {'Ga': 0.8345574460809416})
        self.assertAlmostEqual(get_entropy_from_frequencies({'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}, max_env=2)["Ga"], 0.0)
        self.assertAlmostEqual(get_entropy_from_frequencies({'Ga': ['O:6', 'O:6', 'O:6', 'O:6']}, max_env=2)["Ga"], 1.0)
        self.assertAlmostEqual(
            get_entropy_from_frequencies({'Ga': ['O:6', 'T:4', 'T:4', 'A:2', 'O:6', 'O:6', 'L:2', 'O:6']}, max_env=4)[
                "Ga"], 1.0 - 1.75 / 2.0)
        self.assertAlmostEqual(get_entropy_from_frequencies(
            {'Ga': ['O:6', 'T:4', 'T:4', 'A:2', 'O:6', 'O:6', 'L:2', 'L:2'],
             'Sn': ['O:6', 'T:4', 'T:4', 'A:2', 'O:6', 'O:6', 'L:2', 'O:6']}, max_env=4)["Sn"], 1.0 - 1.75 / 2.0)


class Get_most_frequent_environment_Test(unittest.TestCase):
    def test_get_most_frequent_environment(self):
        self.assertEqual(get_most_frequent_environment({'Ga': ['O:6', 'O:6', 'T:4', 'T:4']})["Ga"], 0.5)
        self.assertDictEqual(get_most_frequent_environment({'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}), {'Ga': 0.5})
        self.assertDictEqual(
            get_most_frequent_environment({'Ga': ['O:6', 'O:6', 'O:6', 'O:6'], 'Sn': ['O:6', 'O:6', 'O:6', 'O:6']}),
            {'Ga': 1.0, 'Sn': 1.0})


class Get_get_mean_CN_from_frequencies_Test(unittest.TestCase):
    def test_get_mean_CN_from_frequencies(self):
        self.assertAlmostEqual(get_mean_CN_from_frequencies({'Ga': ['O:6', 'O:6', 'T:4', 'T:4']})["Ga"], 5)
        self.assertAlmostEqual(get_mean_CN_from_frequencies({'Ga': ['O:6', 'O:6', 'O:6', 'T:4']})["Ga"], 5.5)


class FrequencyEnvironmentPauling1_Test(unittest.TestCase):
    def setUp(self):
        self.matlist = ["mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-306.json", "mp-886.json", "mp-2605.json"]
        self.lse_dict = {}
        self.Pauling_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_dict[mat] = FrequencyEnvironmentPauling1(self.lse_dict[mat])

    def test_get_details(self):
        self.assertDictEqual(self.Pauling_dict["mp-1788.json"].get_details(),
                             {'As': ['T:4', 'T:4', 'T:4', 'T:4', 'O:6', 'O:6', 'O:6', 'O:6']})
        self.assertDictEqual(self.Pauling_dict["mp-7000.json"].get_details(), {'Si': ['T:4', 'T:4', 'T:4']})
        self.assertDictEqual(self.Pauling_dict["mp-19359.json"].get_details(), {'Na': ['O:6'], 'Fe': ['O:6']})
        self.assertDictEqual(self.Pauling_dict["mp-306.json"].get_details(),
                             {'B': ['TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3']})
        self.assertDictEqual(self.Pauling_dict["mp-886.json"].get_details(), {'Ga': ['O:6', 'O:6', 'T:4', 'T:4']})
        self.assertDictEqual(self.Pauling_dict["mp-2605.json"].get_details(), {'Ca': ['O:6']})


class Pauling0_Test(unittest.TestCase):
    def setUp(self):
        self.matlist = ["mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-306.json", "mp-886.json", "mp-2605.json"]
        self.lse_dict = {}
        self.Pauling_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_dict[mat] = Pauling0(self.lse_dict[mat])

    def test_one(self):
        self.assertDictEqual(dict(self.Pauling_dict["mp-1788.json"].get_cations_in_structure()), {'As': 8})
        self.assertDictEqual(dict(self.Pauling_dict["mp-7000.json"].get_cations_in_structure()), {'Si': 3})
        self.assertDictEqual(dict(self.Pauling_dict["mp-19359.json"].get_cations_in_structure()), {'Na': 1, 'Fe': 1})
        self.assertDictEqual(dict(self.Pauling_dict["mp-306.json"].get_cations_in_structure()), {'B': 6})
        self.assertDictEqual(dict(self.Pauling_dict["mp-886.json"].get_cations_in_structure()), {'Ga': 4})
        self.assertDictEqual(dict(self.Pauling_dict["mp-2605.json"].get_cations_in_structure()), {'Ca': 1})


class PaulingRule1_Test(unittest.TestCase):
    """Tests class PaulingRule1 """

    def setUp(self):
        self.matlist = ["mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-306.json", "mp-886.json", "mp-2605.json"]
        self.lse_dict = {}
        self.Pauling_dict = {}
        self.Pauling_dict2 = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_dict[mat] = Pauling1(self.lse_dict[mat])
            self.Pauling_dict2[mat] = Pauling1(self.lse_dict[mat], onlylowerlimit=True)

    def test_get_details(self):
        self.assertDictEqual(self.Pauling_dict["mp-1788.json"].get_details(),
                             {'cat_dependency': {'As': [4, 4]}, 'cat_val_dependency': {'As': {5: [4, 4]}},
                              'Env_fulfilled': 4, 'Env_notfulfilled': 4, 'Env_out_of_list': 0, 'Cat_out_of_list': 0})
        self.assertDictEqual(self.Pauling_dict["mp-7000.json"].get_details(),
                             {'cat_dependency': {'Si': [3, 0]}, 'cat_val_dependency': {'Si': {4: [3, 0]}},
                              'Env_fulfilled': 3, 'Env_notfulfilled': 0, 'Env_out_of_list': 0, 'Cat_out_of_list': 0})
        self.assertDictEqual(self.Pauling_dict["mp-19359.json"].get_details(),
                             {'cat_dependency': {'Na': [1, 0]}, 'cat_val_dependency': {'Na': {1: [1, 0]}},
                              'Env_fulfilled': 1, 'Env_notfulfilled': 0, 'Env_out_of_list': 0, 'Cat_out_of_list': 1})
        self.assertDictEqual(self.Pauling_dict["mp-306.json"].get_details(),
                             {'cat_dependency': {}, 'cat_val_dependency': {}, 'Env_fulfilled': 0, 'Env_notfulfilled': 0,
                              'Env_out_of_list': 6, 'Cat_out_of_list': 0})
        self.assertDictEqual(self.Pauling_dict["mp-886.json"].get_details(),
                             {'cat_dependency': {'Ga': [2, 2]}, 'cat_val_dependency': {'Ga': {3: [2, 2]}},
                              'Env_fulfilled': 2, 'Env_notfulfilled': 2, 'Env_out_of_list': 0, 'Cat_out_of_list': 0})
        self.assertDictEqual(self.Pauling_dict2["mp-886.json"].get_details(),
                             {'cat_dependency': {'Ga': [4, 0]}, 'cat_val_dependency': {'Ga': {3: [4, 0]}},
                              'Env_fulfilled': 4, 'Env_notfulfilled': 0, 'Env_out_of_list': 0, 'Cat_out_of_list': 0})

    def test_is_fulfilled(self):
        self.assertFalse(self.Pauling_dict["mp-1788.json"].is_fulfilled())
        self.assertTrue(self.Pauling_dict["mp-7000.json"].is_fulfilled())
        self.assertFalse(self.Pauling_dict["mp-2605.json"].is_fulfilled())
        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_dict["mp-306.json"].is_fulfilled()

    def test_predict_env_pauling(self):
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.22), ['does not exist'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.25), ['T:4'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.413), ['T:4'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.414), ['O:6'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.591), ['O:6'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.592), ['FO:7'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.644), ['FO:7'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.645), ['SA:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.731), ['SA:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.732), ['TT_1:9', 'C:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(0.99), ['TT_1:9', 'C:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_window(1.0), ['C:12'])

        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.22), [])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.25), ['T:4'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.413), ['T:4'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.414), ['T:4', 'O:6'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.591), ['T:4', 'O:6'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.592),
                         ['T:4', 'O:6', 'FO:7'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.644),
                         ['T:4', 'O:6', 'FO:7'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.645),
                         ['T:4', 'O:6', 'FO:7', 'SA:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.731),
                         ['T:4', 'O:6', 'FO:7', 'SA:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.732),
                         ['T:4', 'O:6', 'FO:7', 'SA:8', 'TT_1:9', 'C:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(0.99),
                         ['T:4', 'O:6', 'FO:7', 'SA:8', 'TT_1:9', 'C:8'])
        self.assertEqual(self.Pauling_dict["mp-306.json"]._predict_env_pauling_lowerlimit(1.0),
                         ['T:4', 'O:6', 'FO:7', 'SA:8', 'TT_1:9', 'C:8', 'C:12'])


class PaulingRule2_Test(unittest.TestCase):
    """Tests class PaulingRule2 """

    def setUp(self):
        self.matlist = ["mp-7000.json", "mp-19418.json"]
        self.lse_dict = {}
        self.Pauling_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_dict[mat] = Pauling2(self.lse_dict[mat])

    def test_get_details(self):
        self.assertDictEqual(self.Pauling_dict["mp-7000.json"].get_details(),
                             {'bvs_for_each_anion': [2.0, 2.0, 2.0, 2.0, 2.0, 2.0],
                              'cations_around_anion': [['Si', 'Si'], ['Si', 'Si'], ['Si', 'Si'], ['Si', 'Si'],
                                                       ['Si', 'Si'], ['Si', 'Si']],
                              'elementwise_fulfillment': {'Si': [12, 0]},
                              'cations_in_structure': Counter({'Si': 3})
                              })
        self.assertDictEqual(self.Pauling_dict["mp-19418.json"].get_details(),
                             {'bvs_for_each_anion': [2.25, 2.25, 2.25, 2.25, 1.75, 1.75, 1.75, 1.75],
                              'cations_around_anion': [['V', 'Cr', 'Cr'], ['V', 'Cr', 'Cr'], ['V', 'Cr', 'Cr'],
                                                       ['V', 'Cr', 'Cr'], ['V', 'Cr'], ['V', 'Cr'], ['V', 'Cr'],
                                                       ['V', 'Cr']],
                              'elementwise_fulfillment': {'Cr': [0, 12], 'V': [0, 8]},
                              'cations_in_structure': Counter({'V': 2, 'Cr': 2})
                              })

    def test_is_fulfilled(self):
        self.assertTrue(self.Pauling_dict["mp-7000.json"].is_fulfilled())
        self.assertFalse(self.Pauling_dict["mp-19418.json"].is_fulfilled())
        self.assertTrue(self.Pauling_dict["mp-19418.json"].is_fulfilled(tolerance=1))



class PaulingRule2_optimized_environments_Test(unittest.TestCase):
    """Tests class PaulingRule2_optimized_environments """

    # TODO: finish the tests and check the test coverage
    def setUp(self):
        #TODO: initilaize normal pauling as well
        self.matlist = ["mp-7000.json", "mp-566090.json","mp-2605.json","mp-4056.json"]
        self.lse_dict = {}
        self.Pauling_dict = {}
        self.Pauling_dict_perc1={}
        self.Pauling_dict_perc02={}
        self.Pauling_dict_realtwo={}

        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_dict[mat] = Pauling2_optimized_environments(self.lse_dict[mat], perc=0.3)
            self.Pauling_dict_perc02[mat] = Pauling2_optimized_environments(self.lse_dict[mat], perc=0.2)
            self.Pauling_dict_perc1[mat] = Pauling2_optimized_environments(self.lse_dict[mat], perc=1.0)
            self.Pauling_dict_realtwo[mat] = Pauling2(self.lse_dict[mat])

    def test_get_symmetry_equivalent_indices(self):
        structure1 = self.lse_dict["mp-7000.json"].structure
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._get_symmetry_eq_indices(structure=structure1),
                             [[0, 1, 2, 3, 4, 5], [6, 7, 8]])
        structure2 = self.lse_dict["mp-566090.json"].structure
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._get_symmetry_eq_indices(structure=structure2),
                             [[0, 1], [2, 7], [3, 5], [4, 6], [8, 10], [9, 12], [11, 15], [13, 14], [16, 29], [17, 39],
                              [18, 23], [19, 34], [20, 25], [21, 26], [22, 27], [24, 36], [28, 37], [30, 35], [31, 33],
                              [32, 38]])

    def test_get_env_indices(self):
        self.assertListEqual(
            self.Pauling_dict["mp-7000.json"]._get_env_indices(self.lse_dict["mp-7000.json"].coordination_environments,
                                                               self.lse_dict["mp-7000.json"].valences, perc=0.0,
                                                               equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]]),
            [[None], [None], [None], [None], [None], [None], [0], [0], [0]])
        self.assertListEqual(
            self.Pauling_dict["mp-7000.json"]._get_env_indices(self.lse_dict["mp-7000.json"].coordination_environments,
                                                               self.lse_dict["mp-7000.json"].valences, perc=0.3,
                                                               equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]]),
            [[None], [None], [None], [None], [None], [None], [0], [0], [0]])
        self.assertListEqual(
            self.Pauling_dict["mp-7000.json"]._get_env_indices(self.lse_dict["mp-7000.json"].coordination_environments,
                                                               self.lse_dict["mp-7000.json"].valences, perc=1.0,
                                                               equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]]),
            [[None], [None], [None], [None], [None], [None], [0], [0], [0]])

        self.assertListEqual(
            self.Pauling_dict["mp-566090.json"]._get_env_indices(
                coordination_environments=self.lse_dict["mp-566090.json"].coordination_environments,
                valences=self.lse_dict["mp-566090.json"].valences, perc=0.8,
                equivalent_indices=[[0, 1], [2, 7], [3, 5], [4, 6],
                                    [8, 10], [9, 12], [11, 15],
                                    [13, 14], [16, 29], [17, 39],
                                    [18, 23], [19, 34], [20, 25],
                                    [21, 26], [22, 27], [24, 36],
                                    [28, 37], [30, 35], [31, 33],
                                    [32, 38]]),
            [[0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [None], [None], [None],
             [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None],
             [None], [None], [None], [None], [None], [None], [None], [None]]
        )

        self.assertListEqual(
            self.Pauling_dict["mp-566090.json"]._get_env_indices(
                coordination_environments=self.lse_dict["mp-566090.json"].coordination_environments,
                valences=self.lse_dict["mp-566090.json"].valences, perc=0.4,
                equivalent_indices=[[0, 1], [2, 7], [3, 5], [4, 6],
                                    [8, 10], [9, 12], [11, 15],
                                    [13, 14], [16, 29], [17, 39],
                                    [18, 23], [19, 34], [20, 25],
                                    [21, 26], [22, 27], [24, 36],
                                    [28, 37], [30, 35], [31, 33],
                                    [32, 38]]),
            [[0], [0], [0], [0], [0], [0], [0], [0], [0, 1], [0], [0, 1], [0], [0], [0], [0], [0], [None], [None],
             [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None],
             [None], [None], [None], [None], [None], [None], [None], [None], [None]]
        )

    def test_get_combinations_cations(self):
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._get_combinations_cations(
            equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]],
            env_indices=[[None], [None], [None], [None], [None], [None], [0], [0], [0]],
            valences=[-2, -2, -2, -2, -2, -2, 4, 4, 4]), [[0]])
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._get_combinations_cations(
            equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]],
            env_indices=[[None], [None], [None], [None], [None], [None], [0, 1], [0, 1], [0, 1]],
            valences=[-2, -2, -2, -2, -2, -2, 4, 4, 4]), [[0], [1]])
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._get_combinations_cations(
            equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]],
            env_indices=[[None], [None], [None], [None], [None], [None], [0, 2], [0, 2], [0, 2]],
            valences=[-2, -2, -2, -2, -2, -2, 4, 4, 4]), [[0], [2]])
        self.assertListEqual(self.Pauling_dict["mp-566090.json"]._get_combinations_cations(
            equivalent_indices=[[0, 1], [2, 7], [3, 5], [4, 6],
                                [8, 10], [9, 12], [11, 15],
                                [13, 14], [16, 29], [17, 39],
                                [18, 23], [19, 34], [20, 25],
                                [21, 26], [22, 27], [24, 36],
                                [28, 37], [30, 35], [31, 33],
                                [32, 38]],
            env_indices=[[0], [0], [0], [0], [0], [0], [0], [0], [0, 1], [0], [0, 1], [0], [0], [0], [0], [0], [None],
                         [None],
                         [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None], [None],
                         [None],
                         [None], [None], [None], [None], [None], [None], [None], [None], [None]],
            valences=[1, 1, 5, 5, 5, 5, 5, 5, 2, 2, 2, 2, 2, 2, 2, 2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2,
                      -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2, -2]),
            [[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0]])

    def test_get_combinations_all_ions(self):
        self.Pauling_dict["mp-7000.json"]._get_combinations_all_ions(combination_array_only_cations=[[0]], equivalent_indices=[[0, 1, 2, 3, 4, 5], [6, 7, 8]], valences=[-2, -2, -2, -2, -2, -2, 4, 4, 4])

    def test_list_combi_recursive(self):
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._list_combi_recursive([[0]],[0,1]),[[0,0],[0,1]])
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._list_combi_recursive([[0,1]],[0,1]),[[0,1,0],[0,1,1]])
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._list_combi_recursive([[0],[1]],[0,1]),[[0,0],[0,1],[1,0],[1,1]])

    def test_get_anion_bvs(self):
        self.assertListEqual(self.Pauling_dict["mp-7000.json"]._get_anions_bvs(),[2.0, 2.0, 2.0, 2.0, 2.0, 2.0])
        self.assertListEqual(self.Pauling_dict_perc1["mp-566090.json"]._get_anions_bvs(),self.Pauling_dict_realtwo["mp-566090.json"]._get_anions_bvs())

    def test_get_details(self):
        self.assertDictEqual(self.Pauling_dict["mp-7000.json"].get_details(10e-2),self.Pauling_dict_realtwo["mp-7000.json"].get_details(10e-2))
        self.assertDictEqual(self.Pauling_dict_perc1["mp-566090.json"].get_details(10e-2),self.Pauling_dict_realtwo["mp-566090.json"].get_details(10e-2))
        self.assertDictEqual(self.Pauling_dict_perc1["mp-4056.json"].get_details(),self.Pauling_dict_perc02["mp-4056.json"].get_details())


    def test_is_fulfilled(self):
        self.assertTrue(self.Pauling_dict["mp-7000.json"].is_fulfilled())
        self.assertFalse(self.Pauling_dict_perc1["mp-566090.json"].is_fulfilled())
        self.assertFalse(self.Pauling_dict["mp-566090.json"].is_fulfilled())



class PaulingConnection_Test(unittest.TestCase):
    """Tests class PaulingConnection """

    def setUp(self):
        self.paulingconnetion = PaulingConnection(DISTANCE=8.0)
        self.matlist = ["mp-7000.json"]
        self.lse_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse

    def test_is_cationic_site(self):
        valences = [1, 2, 3, -2, -2]
        self.assertTrue(self.paulingconnetion._is_cationic_site(0, valences))
        self.assertFalse(self.paulingconnetion._is_cationic_site(3, valences))

    def test_get_oxygen_neighbors(self):
        self.assertEqual(len(self.paulingconnetion._get_oxygen_neighbors(self.lse_dict["mp-7000.json"],
                                                                         site=self.lse_dict["mp-7000.json"].structure[
                                                                             6], r=5, CN=4, ceindex=0)), 4)

    def test_get_cation_neighbors(self):
        self.assertEqual(len(self.paulingconnetion._get_cation_neighbors(self.lse_dict["mp-7000.json"].structure,
                                                                         self.lse_dict["mp-7000.json"].structure[6], 4,
                                                                         self.lse_dict["mp-7000.json"].valences)), 4)

    def test_site_index(self):
        self.assertEqual(self.paulingconnetion._get_site_index(self.lse_dict["mp-7000.json"].structure[6],
                                                               self.lse_dict["mp-7000.json"].structure), 6)


class Pauling3and4_Test(unittest.TestCase):
    """class to test Pauling3and4"""

    def setUp(self):
        self.matlist = ["mp-7000.json", "mp-19418.json"]
        self.lse_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
        self.pauling3and4_dist1 = Pauling3and4()
        self.pauling3and4_dist1.newsetup(self.lse_dict["mp-7000.json"], save_to_file=False, filename=False,
                                         distance=4.5)
        self.pauling3and4_dist2 = Pauling3and4()
        self.pauling3and4_dist2.newsetup(self.lse_dict["mp-7000.json"], save_to_file=False, filename=False)
        self.pauling3and4_dict3 = Pauling3and4()
        self.pauling3and4_dict3.newsetup(self.lse_dict["mp-19418.json"], save_to_file=False, filename=False, distance=4)

    def test_get_connections(self):
        self.assertDictEqual(self.pauling3and4_dist1.get_connections(), {'no': 6, 'corner': 6, 'edge': 0, 'face': 0})
        self.assertDictEqual(self.pauling3and4_dist1.get_connections(maximumCN=4),
                             {'no': 6, 'corner': 6, 'edge': 0, 'face': 0})
        self.assertDictEqual(self.pauling3and4_dist1.get_connections(maximumCN=3),
                             {'corner': 0, 'edge': 0, 'face': 0, 'no': 0})
        self.assertDictEqual(self.pauling3and4_dist2.get_connections(), {'no': 84, 'corner': 6, 'edge': 0, 'face': 0})
        self.assertDictEqual(self.pauling3and4_dict3.get_connections(), {'no': 2, 'corner': 12, 'edge': 2, 'face': 0})
        self.assertDictEqual(self.pauling3and4_dict3.get_connections(maximumCN=20),
                             {'no': 2, 'corner': 12, 'edge': 2, 'face': 0})

    def test_write_read_file(self):
        outfile_path = tempfile.mkstemp(suffix='.json')[1]
        self.pauling3and4 = Pauling3and4()
        self.pauling3and4.newsetup(self.lse_dict["mp-19418.json"], save_to_file=True, filename=outfile_path,
                                   foldername='')
        self.pauling3and4.newsetup(self.lse_dict["mp-19418.json"], save_to_file=True, filename="mp-19418.json",
                                   foldername="tmp_folder")

        self.assertTrue(os.path.isdir("tmp_folder"))

        self.pauling3and4_new = Pauling3and4()
        self.pauling3and4_new.from_file(outfile_path, '')
        self.assertEqual(self.pauling3and4_new.get_connections(), self.pauling3and4.get_connections())

    def tearDown(self):
        if os.path.isdir("tmp_folder"):
            os.remove(os.path.join("tmp_folder", "mp-19418.json"))
            os.rmdir("tmp_folder")


class Pauling3_Test(unittest.TestCase):
    """class to test Pauling3"""

    def setUp(self):
        self.matlist = ["mp-7000.json", "mp-5986.json", "mp-19418.json"]
        self.lse_dict = {}
        self.Pauling_List = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_List[mat] = Pauling3()
            self.Pauling_List[mat].newsetup(lse=lse, save_to_file=False)
        self.Pauling_distance = Pauling3()
        self.Pauling_distance.newsetup(self.lse_dict["mp-19418.json"], save_to_file=False, distance=4)

    def test_is_fulfilled(self):
        self.assertTrue(self.Pauling_List["mp-7000.json"].is_fulfilled())
        self.assertFalse(self.Pauling_List["mp-5986.json"].is_fulfilled())
        self.assertTrue(self.Pauling_List["mp-5986.json"].is_fulfilled(maximumCN=8))
        self.assertFalse(self.Pauling_List["mp-5986.json"].is_fulfilled(maximumCN=12))
        self.assertTrue(self.Pauling_List["mp-5986.json"].is_fulfilled(maximumCN=11))
        self.assertTrue(self.Pauling_List["mp-19418.json"].is_fulfilled(maximumCN=20))
        self.assertTrue(self.Pauling_List["mp-19418.json"].is_fulfilled())
        self.assertTrue(self.Pauling_distance.is_fulfilled(maximumCN=4))
        self.assertTrue(self.Pauling_List["mp-19418.json"].is_fulfilled())

    def test_get_details(self):
        self.assertDictEqual(self.Pauling_List["mp-7000.json"].get_details(maximumCN=None),
                             {'no': 84, 'corner': 6, 'edge': 0, 'face': 0,
                              'species': {'Si': {4: {'no': 168, 'corner': 12, 'edge': 0, 'face': 0}}}})
        self.assertDictEqual(self.Pauling_List["mp-7000.json"].get_details(maximumCN=4),
                             {'no': 84, 'corner': 6, 'edge': 0, 'face': 0,
                              'species': {'Si': {4: {'no': 168, 'corner': 12, 'edge': 0, 'face': 0}}}})
        self.assertDictEqual(self.Pauling_List["mp-7000.json"].get_details(maximumCN=3),
                             {'no': 0, 'corner': 0, 'edge': 0, 'face': 0, 'species': {}})
        self.assertDictEqual(self.Pauling_List["mp-5986.json"].get_details(maximumCN=11),
                             {'no': 10, 'corner': 3, 'edge': 0, 'face': 0,
                              'species': {'Ti': {4: {'no': 20, 'corner': 6, 'edge': 0, 'face': 0}}}})
        self.assertDictEqual(self.Pauling_List["mp-5986.json"].get_details(maximumCN=12),
                             {'no': 38, 'corner': 9, 'edge': 0, 'face': 11,
                              'species': {'Ti': {4: {'no': 44, 'corner': 6, 'edge': 0, 'face': 8}},
                                          'Ba': {2: {'no': 32, 'corner': 12, 'edge': 0, 'face': 14}}}})
        self.assertDictEqual(self.Pauling_List["mp-5986.json"].get_details(),
                             {'no': 38, 'corner': 9, 'edge': 0, 'face': 11,
                              'species': {'Ti': {4: {'no': 44, 'corner': 6, 'edge': 0, 'face': 8}},
                                          'Ba': {2: {'no': 32, 'corner': 12, 'edge': 0, 'face': 14}}}})

        self.assertDictEqual(self.Pauling_List["mp-19418.json"].get_details(),
                             {'no': 96, 'corner': 12, 'edge': 2, 'face': 0,
                              'species': {'V': {5: {'no': 100, 'corner': 12, 'edge': 0, 'face': 0}},
                                          'Cr': {3: {'no': 92, 'corner': 12, 'edge': 4, 'face': 0}}}})
        self.assertDictEqual(self.Pauling_distance.get_details(maximumCN=4),
                             {'no': 2, 'corner': 0, 'edge': 0, 'face': 0,
                              'species': {'V': {5: {'no': 4, 'corner': 0, 'edge': 0, 'face': 0}}}})
        self.assertDictEqual(self.Pauling_List["mp-19418.json"].get_details(maximumCN=10),
                             {'no': 96, 'corner': 12, 'edge': 2, 'face': 0,
                              'species': {'V': {5: {'no': 100, 'corner': 12, 'edge': 0, 'face': 0}},
                                          'Cr': {3: {'no': 92, 'corner': 12, 'edge': 4, 'face': 0}}}})


class Pauling4_Test(unittest.TestCase):
    """class to test Pauling4"""

    def setUp(self):
        # test it for further candidates -> differing in CN, differing in val
        self.matlist = ["mp-7000.json", "mp-19359.json", "mp-1788.json", "mp-5986.json", "mp-19418.json"]
        self.lse_dict = {}
        self.Pauling_List = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_List[mat] = Pauling4()
            self.Pauling_List[mat].newsetup(lse=lse, save_to_file=False)

    def test_is_fulfilled(self):
        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_List["mp-7000.json"].is_fulfilled()
        self.assertTrue(self.Pauling_List["mp-1788.json"].is_fulfilled())
        self.assertFalse(self.Pauling_List["mp-19359.json"].is_fulfilled())
        self.assertFalse(self.Pauling_List["mp-5986.json"].is_fulfilled())

    def test_get_details(self):
        # TODO: make sure this really works for every possible case

        firstdict = self.Pauling_List["mp-7000.json"]._postevaluation4thrule()
        seconddict = self.Pauling_List["mp-7000.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(4, 4, 4,
                                                                                                               4)
        # totally symmetric
        self.assertDictEqual(firstdict["val1:4"]["val2:4"]["CN1:4"]["CN2:4"], seconddict)
        self.assertEqual(firstdict["elementwise"]["Si"]["val1:4"]["val2:4"]["CN1:4"]["CN2:4"]["no"],
                         seconddict["no"] * 2)

        thirddict = self.Pauling_List["mp-19359.json"]._postevaluation4thrule()

        # results for 6 6 13 and 6 6 3 1 should be the same
        self.assertDictEqual(
            self.Pauling_List["mp-19359.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(6, 6, 1, 3),
            self.Pauling_List["mp-19359.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(6, 6, 3, 1))
        self.assertDictEqual(thirddict["val1:1"]["val2:3"]["CN1:6"]["CN2:6"],
                             thirddict["val1:3"]["val2:1"]["CN1:6"]["CN2:6"])

        # two methods should result in the same results
        self.assertDictEqual(thirddict["val1:1"]["val2:3"]["CN1:6"]["CN2:6"],
                             self.Pauling_List["mp-19359.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 1, 3))

        self.assertDictEqual(thirddict["val1:1"]["val2:1"]["CN1:6"]["CN2:6"],
                             self.Pauling_List["mp-19359.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 1, 1))
        self.assertDictEqual(thirddict["val1:3"]["val2:3"]["CN1:6"]["CN2:6"],
                             self.Pauling_List["mp-19359.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 3, 3))

        fourthdict = self.Pauling_List["mp-1788.json"]._postevaluation4thrule()

        # symmetry should be there!
        self.assertDictEqual(fourthdict["val1:5"]["val2:5"]["CN1:4"]["CN2:6"],
                             fourthdict["val1:5"]["val2:5"]["CN1:6"]["CN2:4"])

        self.assertDictEqual(self.Pauling_List["mp-1788.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
            4, 6, 5, 5),
            self.Pauling_List["mp-1788.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                6, 4, 5, 5))
        self.assertDictEqual(fourthdict["val1:5"]["val2:5"]["CN1:4"]["CN2:4"],
                             self.Pauling_List["mp-1788.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 4, 4, 5, 5))
        self.assertDictEqual(fourthdict["val1:5"]["val2:5"]["CN1:4"]["CN2:6"],
                             self.Pauling_List["mp-1788.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 4, 6, 5, 5))

        self.assertDictEqual(fourthdict["val1:5"]["val2:5"]["CN1:6"]["CN2:6"],
                             self.Pauling_List["mp-1788.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 5, 5))
        self.assertDictEqual(fourthdict["val1:5"]["val2:5"]["CN1:4"]["CN2:6"],
                             self.Pauling_List["mp-1788.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 4, 6, 5, 5))

        self.assertDictEqual(self.Pauling_List["mp-1788.json"].get_details(), {'val1:5': {'val2:5': {
            'CN1:4': {'CN2:4': {'no': 44, 'corner': 0, 'edge': 0, 'face': 0},
                      'CN2:6': {'no': 96, 'corner': 16, 'edge': 0, 'face': 0}},
            'CN1:6': {'CN2:4': {'no': 96, 'corner': 16, 'edge': 0, 'face': 0},
                      'CN2:6': {'no': 48, 'corner': 4, 'edge': 0, 'face': 0}}}}, 'elementwise': {'As': {'val1:5': {
            'val2:5': {'CN1:4': {'CN2:4': {'no': 88, 'corner': 0, 'edge': 0, 'face': 0},
                                 'CN2:6': {'no': 192, 'corner': 32, 'edge': 0, 'face': 0}},
                       'CN1:6': {'CN2:4': {'no': 192, 'corner': 32, 'edge': 0, 'face': 0},
                                 'CN2:6': {'no': 96, 'corner': 8, 'edge': 0, 'face': 0}}}}}},
            'maxval': 5, 'minCN': 4})

        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_List["mp-7000.json"].get_details()

        Dict_new = self.Pauling_List["mp-5986.json"].get_details()

        # should not be symmetric
        self.assertDictEqual(Dict_new['val1:4']['val2:2']["CN1:6"]["CN2:12"],
                             self.Pauling_List["mp-5986.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 12, 4, 2))
        self.assertDictEqual(Dict_new['val1:2']['val2:4']["CN1:6"]["CN2:12"],
                             self.Pauling_List["mp-5986.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 12, 2, 4))

        self.assertDictEqual(Dict_new['val1:4']['val2:2']["CN1:12"]["CN2:6"],
                             self.Pauling_List["mp-5986.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 12, 6, 4, 2))

        self.assertDictEqual(Dict_new['val1:4']['val2:4']["CN1:6"]["CN2:6"],
                             self.Pauling_List["mp-5986.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 4, 4))

        self.assertDictEqual(Dict_new['val1:2']['val2:2']["CN1:12"]["CN2:12"],
                             self.Pauling_List["mp-5986.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 12, 12, 2, 2))

        self.assertDictEqual(Dict_new['val1:4']['val2:4']["CN1:6"]["CN2:6"],
                             self.Pauling_List["mp-5986.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 4, 4))
        self.assertDictEqual(Dict_new['val1:4']['val2:2']['CN1:6']['CN2:12'],
                             Dict_new['val1:2']['val2:4']['CN1:12']['CN2:6'])

        Dict2_new = self.Pauling_List["mp-19418.json"].get_details()
        # should not be symmetric:
        self.assertDictEqual(Dict2_new['val1:5']['val2:3']['CN1:4']['CN2:6'],
                             self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 4, 6, 5, 3))
        self.assertDictEqual(self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
            4, 6, 5, 3), {'no': 52, 'corner': 12, 'edge': 0, 'face': 0})
        self.assertDictEqual(self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
            4, 6, 3, 5), {'no': 0, 'corner': 0, 'edge': 0, 'face': 0})
        self.assertDictEqual(Dict2_new['val1:3']['val2:5']['CN1:6']['CN2:4'],
                             self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 4, 3, 5))
        self.assertDictEqual(Dict2_new['val1:3']['val2:5']['CN1:4']['CN2:6'],
                             self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 4, 5, 3))

        self.assertDictEqual(Dict2_new['val1:5']['val2:5']['CN1:4']['CN2:4'],
                             self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 4, 4, 5, 5))

        self.assertDictEqual(Dict2_new['val1:3']['val2:3']['CN1:6']['CN2:6'],
                             self.Pauling_List["mp-19418.json"]._postevaluation4thruleperpolyhedron_only_withoutproduct(
                                 6, 6, 3, 3))


class Pauling5_Test(unittest.TestCase):
    """class to test Pauling5"""

    def setUp(self):
        self.matlist = ["mp-7000.json", "mp-19359.json", "mp-1788.json", "mp-12236.json", "mp-5986.json"]
        self.lse_dict = {}
        self.Pauling_List = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_List[mat] = Pauling5()
            self.Pauling_List[mat].newsetup(lse=lse, save_to_file=False)

    def test_is_fulfilled(self):
        with self.assertRaises(RuleCannotBeAnalyzedError):
            (self.Pauling_List["mp-5986.json"].is_fulfilled())
        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_List["mp-19359.json"].is_fulfilled()
        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_List["mp-7000.json"].is_fulfilled(leave_out_list=["Si"])
        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_List["mp-12236.json"].is_fulfilled(leave_out_list=["Er", "Ga"])

        self.assertFalse(self.Pauling_List["mp-1788.json"].is_fulfilled())
        self.assertFalse(self.Pauling_List["mp-12236.json"].is_fulfilled())
        # TODO: fix problems here

        self.assertTrue(self.Pauling_List["mp-12236.json"].is_fulfilled(leave_out_list=["Ga"]))

        self.assertFalse(self.Pauling_List["mp-12236.json"].is_fulfilled(leave_out_list=["Er"]))
        self.assertTrue(self.Pauling_List["mp-7000.json"].is_fulfilled())

        self.assertFalse(self.Pauling_List["mp-1788.json"].is_fulfilled(options='env'))
        self.assertFalse(self.Pauling_List["mp-12236.json"].is_fulfilled(options='env'))
        self.assertTrue(self.Pauling_List["mp-12236.json"].is_fulfilled(options='env', leave_out_list=["Ga"]))
        self.assertFalse(self.Pauling_List["mp-12236.json"].is_fulfilled(options='env', leave_out_list=["Er"]))
        self.assertTrue(self.Pauling_List["mp-7000.json"].is_fulfilled(options='env'))

        self.assertFalse(self.Pauling_List["mp-1788.json"].is_fulfilled(options='env+nconnections'))
        self.assertFalse(self.Pauling_List["mp-12236.json"].is_fulfilled(options='env+nconnections'))
        self.assertTrue(
            self.Pauling_List["mp-12236.json"].is_fulfilled(options='env+nconnections', leave_out_list=["Ga"]))
        self.assertFalse(
            self.Pauling_List["mp-12236.json"].is_fulfilled(options='env+nconnections', leave_out_list=["Er"]))
        self.assertTrue(self.Pauling_List["mp-7000.json"].is_fulfilled(options='env+nconnections'))

        with self.assertRaises(ValueError):
            self.Pauling_List["mp-7000.json"].is_fulfilled(options='something')

    def test_get_details(self):
        with self.assertRaises(RuleCannotBeAnalyzedError):
            self.Pauling_List["mp-19359.json"].get_details()
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(),
                             {'Ga': {'not_fulfilled': 1, 'fulfilled': 0}, 'Er': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-7000.json"].get_details(),
                             {'Si': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-1788.json"].get_details(),
                             {'As': {'not_fulfilled': 1, 'fulfilled': 0}})
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(options='env'),
                             {'Ga': {'not_fulfilled': 1, 'fulfilled': 0}, 'Er': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-7000.json"].get_details(options='env'),
                             {'Si': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-1788.json"].get_details(options='env'),
                             {'As': {'not_fulfilled': 1, 'fulfilled': 0}})
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(options='env+nconnections'),
                             {'Ga': {'not_fulfilled': 1, 'fulfilled': 0}, 'Er': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-7000.json"].get_details(options='env+nconnections'),
                             {'Si': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-1788.json"].get_details(options='env+nconnections'),
                             {'As': {'not_fulfilled': 1, 'fulfilled': 0}})
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(leave_out_list=["Ga"]),
                             {'Er': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(leave_out_list=["Er"]),
                             {'Ga': {'not_fulfilled': 1, 'fulfilled': 0}})
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(options='env', leave_out_list=["Ga"]),
                             {'Er': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(self.Pauling_List["mp-12236.json"].get_details(options='env', leave_out_list=["Er"]),
                             {'Ga': {'not_fulfilled': 1, 'fulfilled': 0}})
        self.assertDictEqual(
            self.Pauling_List["mp-12236.json"].get_details(options='env+nconnections', leave_out_list=["Ga"]),
            {'Er': {'not_fulfilled': 0, 'fulfilled': 1}})
        self.assertDictEqual(
            self.Pauling_List["mp-12236.json"].get_details(options='env+nconnections', leave_out_list=["Er"]),
            {'Ga': {'not_fulfilled': 1, 'fulfilled': 0}})

        with self.assertRaises(ValueError):
            self.Pauling_List["mp-7000.json"].get_details(options='something')

    def test_write_read_file(self):
        self.Pauling_here = Pauling5()
        outfile_path = tempfile.mkstemp(suffix='.json')[1]
        self.pauling5 = Pauling5()
        self.pauling5.newsetup(self.lse_dict["mp-7000.json"], save_to_file=True, filename=outfile_path, foldername='')

        self.pauling5_new = Pauling5()
        self.pauling5_new.from_file(outfile_path, '')
        self.assertEqual(self.pauling5_new.get_details(), self.pauling5.get_details())

        self.pauling5.newsetup(self.lse_dict["mp-7000.json"], save_to_file=True, filename="mp-7000.json",
                               foldername='tmp_folder')
        self.assertTrue(os.path.isdir("tmp_folder"))

    def tearDown(self):
        if os.path.isdir("tmp_folder"):
            os.remove(os.path.join("tmp_folder", "mp-7000.json"))
            os.rmdir("tmp_folder")


if __name__ == '__main__':
    unittest.main()

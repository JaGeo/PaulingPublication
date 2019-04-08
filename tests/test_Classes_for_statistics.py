from Classes_for_statistics import OverAllAnalysis, Pauling1Frequency, Pauling1Entropy, Pauling1OverAllAnalysis, \
    Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, \
    AllPaulingOverAllAnalysis
from collections import OrderedDict
import unittest
import os
import tempfile
import numpy as np

#
# class TestOverallAnalysis(unittest.TestCase):
#
#     def setUp(self):
#         self.overallanalysis = OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
#                                                lowest_number_environments_for_plot=50, lower_limit_plot=0.0,
#                                                upper_limit_plot=1.0,
#                                                analyse_structures=True, use_prematching=True)
#
#     def test_get_lse_from_folder(self):
#         pass
#
#     def test_get_precomputed_result(self):
#         pass
#
#     def test_add_dict_cat_dependency(self):
#         # test all options for number_of_elements_to_add
#         dict1 = {"Ga": 1, "Sn": 2}
#         dict2 = {"Ga": 1, "Fe": 1}
#         self.overallanalysis._add_dict_cat_dependency(dict1, dict2, number_of_elements_to_add=1)
#         self.assertDictEqual(dict1, {'Ga': 2, 'Sn': 2, 'Fe': 1})
#
#         dict1 = {"Ga": [1, 0], "Sn": [0, 2]}
#         dict2 = {"Ga": [1, 1], "Fe": [1, 0]}
#         self.overallanalysis._add_dict_cat_dependency(dict1, dict2, number_of_elements_to_add=2)
#         self.assertDictEqual(dict1, {'Ga': [2, 1], 'Sn': [0, 2], 'Fe': [1, 0]})
#
#         dict1 = {"Ga": ["O:6", "O:6"], "Sn": ["T:4", "T:4"]}
#         dict2 = {"Ga": ["O:6", "O:6"], "Fe": ["T:4"]}
#
#         self.overallanalysis._add_dict_cat_dependency(dict1, dict2, number_of_elements_to_add=4)
#         self.assertDictEqual(dict1, {'Ga': ["O:6", "O:6", "O:6", "O:6"], 'Sn': ["T:4", "T:4"], 'Fe': ["T:4"]})
#
#     def test_get_similar_structures(self):
#         self.assertDictEqual(dict((self.overallanalysis._get_similar_structures(["mp-7000", "mp-1788", "mp-553432"],
#                                                                                 start_from_Matching=True,
#                                                                                 save_to_file=False)))[
#                                  "structure_matching"], {'mp-7000': ['mp-7000', 'mp-553432'], 'mp-1788': ['mp-1788']})
#
#         self.assertDictEqual(dict(
#             self.overallanalysis._get_similar_structures(["mp-7000", "mp-1788", "mp-553432"], save_to_file=False)[
#                 "structure_matching"]), {'mp-7000': ['mp-7000', 'mp-553432'], 'mp-1788': ['mp-1788']})
#         self.assertDictEqual(dict(
#             self.overallanalysis._get_similar_structures(["mp-7000", "mp-1788", "mp-7000", "mp-553432"],
#                                                          save_to_file=False)["structure_matching"]),
#             {'mp-7000': ['mp-7000', 'mp-553432'], 'mp-1788': ['mp-1788']})
#         outfile_path = tempfile.mkstemp(suffix='.json')[1]
#         self.overallanalysis._get_similar_structures(["mp-7000", "mp-1788", "mp-7000"], save_to_file=True,
#                                                      path_to_save=outfile_path)
#         self.assertDictEqual(
#             self.overallanalysis._get_similar_structures(["mp-553432"], save_to_file=True, path_to_save=outfile_path,
#                                                          restart_from_matching=True),
#             {'list_mat_id': ['mp-7000', 'mp-553432', 'mp-1788'],
#              'structure_matching': OrderedDict([('mp-7000', ['mp-7000', 'mp-553432']), ('mp-1788', ['mp-1788'])]),
#              'additional_info': {'mp-7000': 'SiO2', 'mp-1788': 'As2O5', 'mp-553432': 'TiO2'}})
#         self.assertDictEqual(
#             self.overallanalysis._get_similar_structures(["mp-7000", "mp-1788", "mp-7000", "mp-553432"],
#                                                          fetch_results_only=True, path_to_save=outfile_path),
#             {'list_mat_id': ['mp-7000', 'mp-553432', 'mp-1788'],
#              'structure_matching': OrderedDict([('mp-7000', ['mp-7000', 'mp-553432']), ('mp-1788', ['mp-1788'])]),
#              'additional_info': {'mp-7000': 'SiO2', 'mp-1788': 'As2O5', 'mp-553432': 'TiO2'}})
#         # auch das restarten sollte getestet werden
#
#     def test_print_to_file_similar_structures(self):
#         # TODO: print this stuff to a file and see if it is the way you wanted it to have
#         pass
#
#     def test_pieplot_connections(self):
#         pass
#
#     def test_visualize_structure_by_id(self):
#         pass
#
#
# class TestPauling1Frequency(unittest.TestCase):
#
#     def setUp(self):
#         self.pauling1_normal = Pauling1Frequency(source='my_own_list', onlybinaries=False,
#                                                  plot_element_dependend_analysis=False,
#                                                  list_of_materials_to_investiage='test_list.json')
#         self.pauling1_binary = Pauling1Frequency(source='my_own_list', onlybinaries=True,
#                                                  plot_element_dependend_analysis=False,
#                                                  list_of_materials_to_investiage='test_list.json')
#
#     def test_run(self):
#         self.pauling1_normal.run(save_result_data=False, start_from_results=False, path_to_save='')
#         self.assertDictEqual(self.pauling1_normal.All_Details,
#                              {'As': ['T:4', 'T:4', 'T:4', 'T:4', 'O:6', 'O:6', 'O:6', 'O:6'],
#                               'Si': ['T:4', 'T:4', 'T:4'], 'Na': ['O:6'], 'Fe': ['O:6'],
#                               'B': ['TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3'], 'Ga': ['O:6', 'O:6', 'T:4', 'T:4'],
#                               'Ca': ['O:6']})
#         self.assertDictEqual(self.pauling1_normal.present_env,
#                              {'As': 8, 'Si': 3, 'Na': 1, 'Fe': 1, 'B': 6, 'Ga': 4, 'Ca': 1})
#         self.assertDictEqual(self.pauling1_normal.Plot_PSE_most_frequent,
#                              {'As': 0.5, 'Si': 1.0, 'Na': 1.0, 'Fe': 1.0, 'B': 1.0, 'Ga': 0.5, 'Ca': 1.0})
#
#         self.pauling1_binary.run(save_result_data=False, start_from_results=False, path_to_save='')
#         self.assertDictEqual(self.pauling1_binary.All_Details,
#                              {'As': ['T:4', 'T:4', 'T:4', 'T:4', 'O:6', 'O:6', 'O:6', 'O:6'],
#                               'Si': ['T:4', 'T:4', 'T:4'],
#                               'B': ['TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3'], 'Ga': ['O:6', 'O:6', 'T:4', 'T:4'],
#                               'Ca': ['O:6']})
#         self.assertDictEqual(self.pauling1_binary.present_env,
#                              {'As': 8, 'Si': 3, 'B': 6, 'Ga': 4, 'Ca': 1})
#
#
# class TestPauling1Entropy(unittest.TestCase):
#
#     def setUp(self):
#         self.pauling1_normal = Pauling1Entropy(source='my_own_list', onlybinaries=False,
#                                                plot_element_dependend_analysis=False,
#                                                list_of_materials_to_investiage='test_list.json')
#
#         self.pauling1_binary = Pauling1Entropy(source='my_own_list', onlybinaries=True,
#                                                plot_element_dependend_analysis=False,
#                                                list_of_materials_to_investiage='test_list.json')
#
#     def test_run(self):
#         self.pauling1_normal.run(save_result_data=False, start_from_results=False, path_to_save='')
#         self.assertDictEqual(self.pauling1_normal.All_Details,
#                              {'As': ['T:4', 'T:4', 'T:4', 'T:4', 'O:6', 'O:6', 'O:6', 'O:6'],
#                               'Si': ['T:4', 'T:4', 'T:4'], 'Na': ['O:6'], 'Fe': ['O:6'],
#                               'B': ['TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3'], 'Ga': ['O:6', 'O:6', 'T:4', 'T:4'],
#                               'Ca': ['O:6']})
#         self.assertDictEqual(self.pauling1_normal.present_env,
#                              {'As': 8, 'Si': 3, 'Na': 1, 'Fe': 1, 'B': 6, 'Ga': 4, 'Ca': 1})
#
#         self.assertDictEqual(self.pauling1_normal.Plot_PSE_entropy,
#                              {'As': 0.8345574460809416, 'Si': 1.0, 'Na': 1.0, 'Fe': 1.0, 'B': 1.0,
#                               'Ga': 0.8345574460809416, 'Ca': 1.0})
#
#         self.pauling1_binary.run(save_result_data=False, start_from_results=False, path_to_save='')
#         self.assertDictEqual(self.pauling1_binary.All_Details,
#                              {'As': ['T:4', 'T:4', 'T:4', 'T:4', 'O:6', 'O:6', 'O:6', 'O:6'],
#                               'Si': ['T:4', 'T:4', 'T:4'],
#                               'B': ['TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3', 'TL:3'], 'Ga': ['O:6', 'O:6', 'T:4', 'T:4'],
#                               'Ca': ['O:6']})
#         self.assertDictEqual(self.pauling1_binary.present_env,
#                              {'As': 8, 'Si': 3, 'B': 6, 'Ga': 4, 'Ca': 1})
#
#         # TODO: test save and read file
#
#
# class TestPauling1OverAllAnalysis(unittest.TestCase):
#
#     def setUp(self):
#         self.pauling1_normal = Pauling1OverAllAnalysis(source='my_own_list', onlybinaries=False,
#                                                        plot_element_dependend_analysis=False,
#                                                        list_of_materials_to_investiage='test_list.json',
#                                                        analyse_structures=True, use_prematching=True)
#
#     def test_run(self):
#         self.pauling1_normal.run(start_from_results=False, save_result_data=False, path_to_save='',
#                                  save_structure_analysis=False)
#         self.assertDictEqual(self.pauling1_normal.Plot_PSE_DICT,
#                              {'As': [4, 4], 'Si': [3, 0], 'Na': [1, 0], 'Ga': [2, 2], 'Ca': [0, 1]})
#         self.assertListEqual(self.pauling1_normal.structures_cannot_be_evaluated, ['mp-19359', 'mp-306'])
#         self.assertListEqual(self.pauling1_normal.structures_exceptions, ['mp-1788', 'mp-886', 'mp-2605'])
#         self.assertListEqual(self.pauling1_normal.structures_fulfillingrule, ['mp-7000'])
#         self.assertDictEqual(self.pauling1_normal.dict_similarstructures_fulfilling,
#                              {'list_mat_id': ['mp-7000'], 'structure_matching': OrderedDict([('mp-7000', ['mp-7000'])]),
#                               'additional_info': {'mp-7000': 'SiO2'}})
#         self.assertDictEqual(self.pauling1_normal.dict_similarstructures_exceptions,
#                              {'list_mat_id': ['mp-1788', 'mp-886', 'mp-2605'], 'structure_matching': OrderedDict(
#                                  [('mp-1788', ['mp-1788']), ('mp-886', ['mp-886']), ('mp-2605', ['mp-2605'])]),
#                               'additional_info': {'mp-1788': 'As2O5', 'mp-886': 'Ga2O3', 'mp-2605': 'CaO'}})
#
#         # TODO: printing to file and to screen
#
#
# class TestPauling2OverAllAnalysis(unittest.TestCase):
#
#     def setUp(self):
#         self.pauling2_normal = Pauling2OverAllAnalysis(source='my_own_list', onlybinaries=False,
#                                                        plot_element_dependend_analysis=False,
#                                                        list_of_materials_to_investiage='test_list.json',
#                                                        analyse_structures=True, use_prematching=True)
#
#     def test_run(self):
#         self.pauling2_normal.run(start_from_results=False, save_result_data=False, path_to_save='',
#                                  show_plot=False,
#                                  save_structure_analysis=False)
#         self.assertListEqual(self.pauling2_normal.structures_fulfillingrule,
#                              ['mp-7000', 'mp-19359', 'mp-306', 'mp-2605'])
#         self.assertListEqual(self.pauling2_normal.structures_exceptions, ["mp-1788", "mp-886"])
#         self.assertListEqual(self.pauling2_normal.structures_cannot_be_evaluated, [])
#         self.assertDictEqual(self.pauling2_normal.dict_similarstructures_exceptions,
#                              {'list_mat_id': ['mp-1788', 'mp-886'],
#                               'structure_matching': OrderedDict([('mp-1788', ['mp-1788']), ('mp-886', ['mp-886'])]),
#                               'additional_info': {'mp-1788': 'As2O5', 'mp-886': 'Ga2O3'}})
#         self.assertDictEqual(self.pauling2_normal.dict_similarstructures_fulfilling,
#                              {'list_mat_id': ['mp-7000', 'mp-19359', 'mp-306', 'mp-2605'],
#                               'structure_matching': OrderedDict(
#                                   [('mp-7000', ['mp-7000']), ('mp-19359', ['mp-19359']), ('mp-306', ['mp-306']),
#                                    ('mp-2605', ['mp-2605'])]),
#                               'additional_info': {'mp-7000': 'SiO2', 'mp-19359': 'NaFeO2', 'mp-306': 'B2O3',
#                                                   'mp-2605': 'CaO'}})
#         self.assertDictEqual(self.pauling2_normal.present_env,
#                              {'As': 8, 'Si': 3, 'Na': 1, 'Fe': 1, 'B': 6, 'Ga': 4, 'Ca': 1})
#         self.assertAlmostEqual(self.pauling2_normal.tot_stddev[0],
#                                np.std(self.pauling2_normal.array_bvs, ddof=1) / np.mean(self.pauling2_normal.array_bvs))
#         # charge neutrality:
#         self.assertEqual(self.pauling2_normal.bs_sum_mean, 2.0)
#
#         for ivalue, value in enumerate(self.pauling2_normal.array_bvs):
#             if value > 2.0:
#                 newvalue = value - self.pauling2_normal.bs_dev[ivalue]
#             else:
#                 newvalue = value + self.pauling2_normal.bs_dev[ivalue]
#             self.assertAlmostEqual(newvalue, 2.0)
#
#         comparedict2 = [0., 0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065,
#                         0.07, 0.075, 0.08, 0.085, 0.09, 0.095, 0.1, 0.105, 0.11, 0.115, 0.12, 0.125, 0.13, 0.135,
#                         0.14, 0.145, 0.15, 0.155, 0.16, 0.165, 0.17]
#         for ivalue, value in enumerate(list(self.pauling2_normal.arraydev_share)):
#             self.assertAlmostEqual(value, comparedict2[ivalue])
#
#         comparedic3 = [0.43181818, 0.45454545, 0.45454545, 0.45454545, 0.45454545, 0.45454545, 0.45454545, 0.45454545,
#                        0.45454545, 0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182,
#                        0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182, 0.81818182,
#                        0.81818182, 0.90909091, 0.90909091, 0.90909091, 0.90909091, 0.90909091, 0.90909091, 0.90909091,
#                        0.90909091, 0.90909091, 1.]
#         for ivalue, value in enumerate(list(self.pauling2_normal.relativefrequency)):
#             self.assertAlmostEqual(value, comparedic3[ivalue])
#
#         # todo: read and write file
#
#     def test_stddev(self):
#         # test case from https://en.wikipedia.org/wiki/Standard_deviation
#         list_here = [727.7, 1086.5, 1091.0, 1361.3, 1490.5, 1956.1]
#         mean = 1285.5
#         result = 420.96 / 1285.5
#         self.assertAlmostEqual(self.pauling2_normal._stddev(list_here, mean=mean)[0], result, 4)
#
#     def test_get_deviation_from_ideal_value(self):
#         list_here = [0.5, 1.0, -2, 3, 0]
#         ideal = 0
#         self.assertListEqual(self.pauling2_normal._get_deviation_from_ideal_value(list_here, ideal),
#                              [0.5, 1.0, 2.0, 3.0, 0.0])
#         ideal = 1.0
#         self.assertListEqual(self.pauling2_normal._get_deviation_from_ideal_value(list_here, ideal),
#                              [0.5, 0.0, 3.0, 2.0, 1.0])
#
#     def test_get_frequency_of_values(self):
#         comparelist = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
#         for iel, el in enumerate(list(self.pauling2_normal._get_frequency_of_values([0.0, 0.5, 0.7], 0.1)[0])):
#             self.assertAlmostEqual(el, comparelist[iel])
#
#         self.assertListEqual(list(self.pauling2_normal._get_frequency_of_values([0.0, 0.5, 0.7], 0.1)[1]),
#                              [1, 1, 1, 1, 1, 2, 2, 3])
#
#         self.assertListEqual(list(self.pauling2_normal._get_frequency_of_values([0.1, 0.5, 0.6], 0.1)[1]),
#                              [0, 1, 1, 1, 1, 2, 3])
#
#     def test_list_to_np_array_and_divide_by_value(self):
#         list_here = [2.1, 2.1, 1.9, 1]
#         divisor = 2
#         self.assertListEqual(list(self.pauling2_normal._list_to_np_array_and_divide_by_value(list_here, divisor)),
#                              [1.05, 1.05, 0.95, 0.5])
#
#
# class TestPauling3OverAllAnalysis(unittest.TestCase):
#     def setUp(self):
#         self.pauling3normal = Pauling3OverAllAnalysis(source='my_own_list', onlybinaries=False,
#                                                       plot_element_dependend_analysis=False,
#                                                       list_of_materials_to_investiage='test_list.json',
#                                                       analyse_structures=True, use_prematching=True)
#
#     def test_run(self):
#         self.pauling3normal.run(show_plot=False, start_from_connections=False, save_connections=False,
#                                 connections_folder='', start_from_results=False, save_result_data=False,
#                                 restart_from_saved_structure_analyisis=False, save_structure_analysis=False,
#                                 path_to_save='')
#         self.assertDictEqual(self.pauling3normal.present_env,
#                              {'As': 8, 'Si': 3, 'Na': 1, 'Fe': 1, 'B': 6, 'Ga': 4, 'Ca': 1})
#         allcorneredge = 0
#         allface = 0
#         for item in self.pauling3normal.Plot_PSE_DICT.values():
#             allcorneredge += item[0]
#             allface += item[1]
#         self.assertEqual(allcorneredge,
#                          (self.pauling3normal.connections['corner'] + self.pauling3normal.connections['edge']) * 2)
#
#         comparisondict = {"corner": 20 + 6 + 6 + 9 + 16 + 3, "edge": 12 + 4 + 6, "face": 0}
#         self.assertDictEqual(comparisondict, self.pauling3normal.connections)
#         self.assertEqual(allface, (self.pauling3normal.connections['face']) * 2)
#
#         # self.assertDictEqual()
#         self.assertListEqual(self.pauling3normal.structures_cannot_be_evaluated, [])
#         self.assertListEqual(self.pauling3normal.structures_exceptions, [])
#         self.assertListEqual(self.pauling3normal.structures_fulfillingrule,
#                              ['mp-1788', 'mp-7000', 'mp-19359', 'mp-306', 'mp-886', 'mp-2605'])
#
#         self.assertDictEqual(self.pauling3normal.dict_similarstructures_fulfilling,
#                              {'list_mat_id': ['mp-1788', 'mp-7000', 'mp-19359', 'mp-306', 'mp-886', 'mp-2605'],
#                               'structure_matching': OrderedDict(
#                                   [('mp-1788', ['mp-1788']), ('mp-7000', ['mp-7000']), ('mp-19359', ['mp-19359']),
#                                    ('mp-306', ['mp-306']), ('mp-886', ['mp-886']), ('mp-2605', ['mp-2605'])]),
#                               'additional_info': {'mp-1788': 'As2O5', 'mp-7000': 'SiO2', 'mp-19359': 'NaFeO2',
#                                                   'mp-306': 'B2O3', 'mp-886': 'Ga2O3', 'mp-2605': 'CaO'}})
#         self.assertDictEqual(self.pauling3normal.dict_similarstructures_exceptions,
#                              {'list_mat_id': [], 'structure_matching': OrderedDict(), 'additional_info': {}})
#
#         # TODO: more tests for reformat etc
#
#
# class TestPauling4OverAllAnalysis(unittest.TestCase):
#     def setUp(self):
#         self.pauling4normal = Pauling4OverAllAnalysis(source='my_own_list', onlybinaries=False,
#                                                       plot_element_dependend_analysis=False,
#                                                       list_of_materials_to_investiage='test_list.json',
#                                                       analyse_structures=True, use_prematching=True)
#
#     def test_run(self):
#         self.pauling4normal.run(show_plot=False, start_from_connections=False, save_connections=False,
#                                 connections_folder='', start_from_results=False, save_result_data=False,
#                                 restart_from_saved_structure_analyisis=False, save_structure_analysis=False,
#                                 path_to_save='')
#
#         self.assertDictEqual(self.pauling4normal.Plot_PSE_DICT,{'As': [88, 0], 'Fe': [36, 6], 'Ga': [76, 4]})
#         self.assertDictEqual(self.pauling4normal.present_env,{'As': 8, 'Na': 1, 'Fe': 1, 'Ga': 4})
#         self.assertListEqual(self.pauling4normal.structures_fulfillingrule,['mp-1788'])
#         self.assertListEqual(self.pauling4normal.structures_exceptions,['mp-19359', 'mp-886'])
#         self.assertListEqual(self.pauling4normal.structures_cannot_be_evaluated,['mp-7000', 'mp-306', 'mp-2605'])
#
#         self.assertDictEqual(self.pauling4normal.Dict_val,{'5': {'5': [284, 36]}, '1': {'1': [18, 3], '3': [60, 24]}, '3': {'1': [60, 24], '3': [210, 37]}})
#
#         self.assertDictEqual(self.pauling4normal.Dict_CN,{'4': {'4': [82, 2], '6': [308, 60]}, '6': {'4': [308, 60], '6': [182, 38]}})
#
#     def test_reformatting(self):
#         self.assertDictEqual(self.pauling4normal._reformat_details_CN({'val1:3': {'val2:3': {'CN1:6': {'CN2:6': {'no': 38, 'corner': 2, 'edge': 4, 'face': 3}, 'CN2:4': {'no': 58, 'corner': 14, 'edge': 2, 'face': 2}}, 'CN1:4': {'CN2:6': {'no': 58, 'corner': 14, 'edge': 1, 'face': 1}, 'CN2:4': {'no': 38, 'corner': 2, 'edge': 1, 'face': 1}}}}}),{"6":{"6": [38,9], "4": [116,34] }, "4": {"6": [116,34], "4": [38,4]}})
#         self.assertDictEqual(self.pauling4normal._reformat_details_val({'val1:3': {'val2:3': {'CN1:6': {'CN2:6': {'no': 38, 'corner': 2, 'edge': 4, 'face': 3}, 'CN2:4': {'no': 58, 'corner': 14, 'edge': 2, 'face': 2}}, 'CN1:4': {'CN2:6': {'no': 58, 'corner': 14, 'edge': 1, 'face': 1}, 'CN2:4': {'no': 38, 'corner': 2, 'edge': 1, 'face': 1}}}}}),{'3': {'3': [192, 47]}})
#         self.assertDictEqual(self.pauling4normal._reformat_details_val({'val1:1': {'val2:1': {'CN1:6': {'CN2:6': {'no': 18, 'corner': 0, 'edge': 3, 'face': 0}}}, 'val2:3': {'CN1:6': {'CN2:6': {'no': 30, 'corner': 6, 'edge': 6, 'face': 0}}}}, 'val1:3': {'val2:1': {'CN1:6': {'CN2:6': {'no': 30, 'corner': 6, 'edge': 6, 'face': 0}}}, 'val2:3': {'CN1:6': {'CN2:6': {'no': 18, 'corner': 0, 'edge': 3, 'face': 0}}}}}),{'1': {'1': [18, 3], '3': [60, 24]}, '3': {'1': [60, 24], '3': [18, 3]}})
#         self.assertDictEqual(self.pauling4normal._reformat_details_elementwise({"maxval":3, "minCN": 4, "elementwise": {'Ga': {'val1:3': {'val2:3': {'CN1:6': {'CN2:6': {'no': 76, 'corner': 0, 'edge': 8, 'face': 0}, 'CN2:4': {'no': 116, 'corner': 28, 'edge': 0, 'face': 0}}, 'CN1:4': {'CN2:6': {'no': 116, 'corner': 28, 'edge': 0, 'face': 0}, 'CN2:4': {'no': 76, 'corner': 4, 'edge': 1, 'face': 1}}}}}}}),{'Ga': [76, 6]})
#
# class TestPauling5OverAllAnalysis(unittest.TestCase):
#     def setUp(self):
#         self.pauling5normal = Pauling5OverAllAnalysis(source='my_own_list', onlybinaries=False,
#                                                       plot_element_dependend_analysis=False,
#                                                       list_of_materials_to_investiage='test_list.json',
#                                                       analyse_structures=True, use_prematching=True)
#
#         self.pauling5entropy=Pauling5OverAllAnalysis(source='my_own_list', onlybinaries=False,
#                                                       plot_element_dependend_analysis=False,
#                                                       list_of_materials_to_investiage='test_list.json',
#                                                       analyse_structures=True, use_prematching=True)
#
#
#     def test_run(self):
#         self.pauling5normal.run(show_plot=False, start_from_connections=False, save_connections=False,
#                                 connections_folder='', start_from_results=False, save_result_data=False,
#                                 restart_from_saved_structure_analyisis=False, save_structure_analysis=False,
#                                 path_to_save='',remove_elements_low_entropy=False)
#         self.assertListEqual(self.pauling5normal.structures_exceptions,['mp-1788', 'mp-886'])
#         self.assertListEqual(self.pauling5normal.structures_fulfillingrule,['mp-7000', 'mp-306'])
#         self.assertListEqual(self.pauling5normal.structures_cannot_be_evaluated,['mp-19359', 'mp-2605'])
#
#         self.assertDictEqual(self.pauling5normal.Plot_PSE_DICT,{'As': [0, 1], 'Si': [1, 0], 'B': [1, 0], 'Ga': [0, 1]})
#         self.assertDictEqual(self.pauling5normal.present_env,{'As': 8, 'Si': 3, 'B': 6, 'Ga': 4})
#
#
#         self.pauling5entropy.run(show_plot=False, start_from_connections=False, save_connections=False,
#                                 connections_folder='', start_from_results=False, save_result_data=False,
#                                 restart_from_saved_structure_analyisis=False, save_structure_analysis=False,
#                                 path_to_save='',remove_elements_low_entropy=True, threshold_remove_elements=0.84)
#
#         self.assertListEqual(self.pauling5entropy.structures_fulfillingrule,[])
#         self.assertListEqual(self.pauling5entropy.structures_exceptions,['mp-1788', 'mp-886'])
#         self.assertListEqual(self.pauling5entropy.structures_cannot_be_evaluated,['mp-7000', 'mp-19359', 'mp-306', 'mp-2605'])
#
#
#         self.assertDictEqual(self.pauling5entropy.Plot_PSE_DICT, {'As': [0, 1], 'Ga': [0, 1]})
#         self.assertDictEqual(self.pauling5entropy.present_env, {'As': 8, 'Ga': 4})
#         self.assertListEqual(self.pauling5entropy.list_to_remove,['Si', 'Na', 'Fe', 'B', 'Ca'])
#
#
#         self.pauling5entropy.run(show_plot=False, start_from_connections=False, save_connections=False,
#                                 connections_folder='', start_from_results=False, save_result_data=False,
#                                 restart_from_saved_structure_analyisis=False, save_structure_analysis=False,
#                                 path_to_save='',remove_elements_low_entropy=True, threshold_remove_elements=0.82)
#
#         self.assertListEqual(self.pauling5entropy.structures_fulfillingrule,[])
#         self.assertListEqual(self.pauling5entropy.structures_exceptions,[])
#         self.assertListEqual(self.pauling5entropy.structures_cannot_be_evaluated,['mp-1788', 'mp-7000', 'mp-19359', 'mp-306','mp-886', 'mp-2605'])
#
#         self.assertDictEqual(self.pauling5entropy.Plot_PSE_DICT, {})
#         self.assertDictEqual(self.pauling5entropy.present_env, {})
#         self.assertListEqual(self.pauling5entropy.list_to_remove, ['As', 'Si', 'Na', 'Fe', 'B', 'Ga', 'Ca'])
#
#     def test_reformat(self):
#         self.assertDictEqual(self.pauling5normal._reformat_details_elementdependency({'As': {'not_fulfilled': 1, 'fulfilled': 0}}),{'As': [0, 1]} )
#         self.assertDictEqual(
#             self.pauling5normal._reformat_details_elementdependency({'As': {'not_fulfilled': 0, 'fulfilled': 1}}),
#             {'As': [1, 0]})
#


class TestAllPaulingOverAllAnalysis(unittest.TestCase):
    #TODO: implmement this testcase and make sure it works for everything and every setting!
    #TODO: implement class in correct manner
    def setUp(self):
        self.paulingall=AllPaulingOverAllAnalysis(source='my_own_list', onlybinaries=False,
                                                       plot_element_dependend_analysis=False,
                                                       list_of_materials_to_investiage='test_list.json',
                                                       analyse_structures=True, use_prematching=True)
        if not os.path.isdir("tmp_folder3"):
            os.mkdir("tmp_folder3")

    def test_run(self):
        # self.paulingall.run(remove_elements_low_entropy=False, start_from_connections=False,
        #     save_connections=False, connections_folder34='AnalysisConnections', connections_folder5='AnalysisConnections_5thRule',
        #     start_from_results=False, save_result_data=False,
        #     restart_from_saved_structure_analyisis=False, save_structure_analysis=False,
        #     path_to_save='', start_material=None, stop_material=None,
        #     threshold_remove_elements=0.95)
        #
        # self.assertListEqual(self.paulingall.structures_cannot_be_evaluated,['mp-7000', 'mp-19359', 'mp-306'])
        # self.assertListEqual(self.paulingall.structures_fulfillingrule,[])
        # self.assertListEqual(self.paulingall.structures_exceptions,['mp-1788', 'mp-886', 'mp-2605'])

        self.paulingall.run(remove_elements_low_entropy=False,start_from_connections=False,
                            save_connections=True,connections_folder34='tmp_folder1',connections_folder5='tmp_folder2',
                            start_from_results=False,save_result_data=False,path_to_save='Results.json', save_structure_analysis=False,
                            start_material=None,stop_material=1
                            )

        self.paulingall.run(remove_elements_low_entropy=False, start_from_connections=True,
                            save_connections=True, connections_folder34='tmp_folder1', connections_folder5='tmp_folder2',
                            start_from_results=False, save_result_data=True, path_to_save='Results.json',
                            save_structure_analysis=False,start_material=None,stop_material=1
                            )


        self.paulingall.run(remove_elements_low_entropy=False, start_from_connections=True,
                           save_connections=True, connections_folder34='tmp_folder1', connections_folder5='tmp_folder2',
                           start_from_results=True, save_result_data=False, path_to_save='Results.json',
                           save_structure_analysis=True,start_material=None,stop_material=1
                           )
        self.assertListEqual(self.paulingall.structures_cannot_be_evaluated,[])
        self.assertListEqual(self.paulingall.structures_fulfillingrule,[])
        self.assertListEqual(self.paulingall.structures_exceptions,['mp-1788'])



        #test reading and writing of all files
        #do that with teardown at the end and clean all files

    def tearDown(self):
        os.remove("Results.json")
        os.remove("Results_structural_exceptions.json")
        os.remove("Results_structural_exceptions_readable.csv")
        os.remove("Results_structural_exceptions_readable.yaml")
        os.remove("Results_structural_fulfilling_readable.yaml")
        os.remove("Results_structural_fulfilling_readable.csv")
        os.remove("Results_structural_fulfilling.json")

        if os.path.isdir("tmp_folder1"):
            os.rmdir("tmp_folder1")

        if os.path.isdir("tmp_folder2"):
            os.rmdir("tmp_folder2")


if __name__ == '__main__':
    unittest.main()

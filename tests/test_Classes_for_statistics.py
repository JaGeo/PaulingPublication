from Classes_for_statistics import OverAllAnalysis, Pauling1Frequency, Pauling1Entropy, Pauling1OverAllAnalysis,Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, AllPaulingOverAllAnalysis
from collections import OrderedDict
import unittest
import tempfile


class TestOverallAnalysis(unittest.TestCase):

    def setUp(self):
        self.overallanalysis= OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 analyse_structures=True, use_prematching=True)


    def test_get_lse_from_folder(self):
        pass

    def test_get_precomputed_result(self):
        pass

    def test_add_dict_cat_dependency(self):
        #test all options for number_of_elements_to_add
        dict1={"Ga":1, "Sn":2}
        dict2={"Ga":1,"Fe":1}
        self.overallanalysis._add_dict_cat_dependency(dict1,dict2,number_of_elements_to_add=1)
        self.assertDictEqual(dict1,{'Ga': 2, 'Sn': 2, 'Fe': 1})

        dict1 = {"Ga": [1,0], "Sn": [0,2]}
        dict2 = {"Ga": [1,1], "Fe": [1,0]}
        self.overallanalysis._add_dict_cat_dependency(dict1,dict2,number_of_elements_to_add=2)
        self.assertDictEqual(dict1,{'Ga': [2,1], 'Sn': [0,2], 'Fe': [1,0]})



        dict1={"Ga":["O:6","O:6"], "Sn":["T:4","T:4"]}
        dict2={"Ga":["O:6","O:6"], "Fe":["T:4"]}

        self.overallanalysis._add_dict_cat_dependency(dict1, dict2, number_of_elements_to_add=4)
        self.assertDictEqual(dict1, {'Ga': ["O:6","O:6","O:6","O:6"], 'Sn': ["T:4","T:4"], 'Fe': ["T:4"]})

    def test_get_similar_structures(self):
        self.assertDictEqual(dict((self.overallanalysis._get_similar_structures(["mp-7000", "mp-1788", "mp-553432"], start_from_Matching=True,save_to_file=False)))["structure_matching"],{'mp-7000': ['mp-7000', 'mp-553432'], 'mp-1788': ['mp-1788']})

        self.assertDictEqual(dict(self.overallanalysis._get_similar_structures(["mp-7000","mp-1788","mp-553432"],save_to_file=False)["structure_matching"]),{'mp-7000': ['mp-7000', 'mp-553432'], 'mp-1788': ['mp-1788']})
        self.assertDictEqual(dict(self.overallanalysis._get_similar_structures(["mp-7000","mp-1788","mp-7000","mp-553432"],save_to_file=False)["structure_matching"]),{'mp-7000': ['mp-7000', 'mp-553432'], 'mp-1788': ['mp-1788']})
        outfile_path = tempfile.mkstemp(suffix='.json')[1]
        self.overallanalysis._get_similar_structures(["mp-7000","mp-1788","mp-7000"],save_to_file=True,path_to_save=outfile_path)
        self.assertDictEqual(self.overallanalysis._get_similar_structures(["mp-553432"],save_to_file=True,path_to_save=outfile_path,restart_from_matching=True),{'list_mat_id': ['mp-7000', 'mp-553432', 'mp-1788'], 'structure_matching': OrderedDict([('mp-7000', ['mp-7000', 'mp-553432']), ('mp-1788', ['mp-1788'])]), 'additional_info': {'mp-7000': 'SiO2', 'mp-1788': 'As2O5', 'mp-553432': 'TiO2'}})
        self.assertDictEqual(self.overallanalysis._get_similar_structures(["mp-7000","mp-1788","mp-7000","mp-553432"],fetch_results_only=True,path_to_save=outfile_path),{'list_mat_id': ['mp-7000', 'mp-553432', 'mp-1788'], 'structure_matching': OrderedDict([('mp-7000', ['mp-7000', 'mp-553432']), ('mp-1788', ['mp-1788'])]), 'additional_info': {'mp-7000': 'SiO2', 'mp-1788': 'As2O5', 'mp-553432': 'TiO2'}})
        #auch das restarten sollte getestet werden


    def test_print_to_file_similar_structures(self):
        #TODO: print this stuff to a file and see if it is the way you wanted it to have
        pass

    def test_pieplot_connections(self):
        pass

    def test_visualize_structure_by_id(self):
        pass


class Pauling1FrequencyTest(unittest.TestCase):
    def setUp(self):
        pass

    def test_(self):
        pass

    #teste klasse ausfuehrlich und auch analysen bis hin zu den plots







if __name__ == '__main__':
    unittest.main()

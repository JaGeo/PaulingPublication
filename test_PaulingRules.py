from PaulingRules import Pauling1, Pauling2, PaulingConnection, Pauling3and4
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import unittest
import json


class PaulingRule1_Test(unittest.TestCase):
    """Tests class PaulingRule1 """

    def setUp(self):
        self.matlist = ["mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-306.json", "mp-886.json"]
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


class PaulingRule2_Test(unittest.TestCase):
    """Tests class PaulingRule2 """

    def setUp(self):
        self.matlist = ["mp-7000.json","mp-19418.json"]
        self.lse_dict = {}
        self.Pauling_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
            self.Pauling_dict[mat] = Pauling2(self.lse_dict[mat])

    def test_get_details(self):
        self.assertDictEqual(self.Pauling_dict["mp-7000.json"].get_details(),{'bvs_for_each_anion': [2.0, 2.0, 2.0, 2.0, 2.0, 2.0], 'cations_around_anion': [['Si', 'Si'], ['Si', 'Si'], ['Si', 'Si'], ['Si', 'Si'], ['Si', 'Si'], ['Si', 'Si']]})
        self.assertDictEqual(self.Pauling_dict["mp-19418.json"].get_details(),{'bvs_for_each_anion': [2.25, 2.25, 2.25, 2.25, 1.75, 1.75, 1.75, 1.75], 'cations_around_anion': [['V', 'Cr', 'Cr'], ['V', 'Cr', 'Cr'], ['V', 'Cr', 'Cr'], ['V', 'Cr', 'Cr'], ['V', 'Cr'], ['V', 'Cr'], ['V', 'Cr'], ['V', 'Cr']]})

    def test_is_fulfilled(self):
        self.assertTrue(self.Pauling_dict["mp-7000.json"].is_fulfilled())
        self.assertFalse(self.Pauling_dict["mp-19418.json"].is_fulfilled())


class PaulingConnection_Test(unittest.TestCase):
    """Tests class PaulingConnection """
    def setUp(self):
        self.paulingconnetion=PaulingConnection(DISTANCE=8.0)
        self.matlist = ["mp-7000.json"]
        self.lse_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse

    def test_is_cationic_site(self):
        valences=[1,2,3,-2,-2]
        self.assertTrue(self.paulingconnetion._is_cationic_site(0, valences))
        self.assertFalse(self.paulingconnetion._is_cationic_site(3, valences))

    def test_get_oxygen_neighbors(self):
        self.assertEqual(len(self.paulingconnetion._get_oxygen_neighbors(self.lse_dict["mp-7000.json"],site=self.lse_dict["mp-7000.json"].structure[6],r=5,CN=4,ceindex=0)),4)

    def test_get_cation_neighbors(self):
        self.assertEqual(len(self.paulingconnetion._get_cation_neighbors(self.lse_dict["mp-7000.json"].structure,self.lse_dict["mp-7000.json"].structure[6],4,self.lse_dict["mp-7000.json"].valences)),4)

    def test_site_index(self):
        self.assertEqual(self.paulingconnetion._get_site_index(self.lse_dict["mp-7000.json"].structure[6],self.lse_dict["mp-7000.json"].structure),6)


class Pauling3and4_Test(unittest.TestCase):
    """class to test Pauling3and4"""
    def setUp(self):

        self.matlist = ["mp-7000.json"]
        self.lse_dict = {}
        for mat in self.matlist:
            with open(mat, "r") as f:
                dict_lse = json.load(f)
            lse = LightStructureEnvironments.from_dict(dict_lse)
            self.lse_dict[mat] = lse
        self.pauling3and4_dist1 = Pauling3and4()
        self.pauling3and4_dist1.newsetup(self.lse_dict["mp-7000.json"],save_to_file=False,filename=False,distance=4.5)
        self.pauling3and4_dist2 = Pauling3and4()
        self.pauling3and4_dist2.newsetup(self.lse_dict["mp-7000.json"],save_to_file=False,filename=False)

    def test(self):
        pass




#
# #"mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-5986.json", "mp-306.json"
# with open("mp-306.json","r") as f:
#     dict_lse=json.load(f
#
#
# lse=LightStructureEnvironments.from_dict(dict_lse)
# print(lse.structure)
#
# pauling1=Pauling1(lse)
#
# #schoenes ausgabeformat überlegen
# print(pauling1.get_details())
#
# print(pauling1.is_fulfilled())
#
# #write a nice unit test - running through nearly all parts of the code
#
#
#
# #brauche lse, um zu testen

if __name__ == '__main__':
    unittest.main()

import json
import math
import os
from collections import Counter
from collections import OrderedDict

import numpy as np
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.core import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


class RuleCannotBeAnalyzedError(Exception):
    """
    An Exception that will be raised if one Rule cannot be analyzed
    """

    def __init__(self, value='The Rule cannot be analyzed'):
        """

        :param value: string that is used in the string representation of the error
        """
        self.value = value

    def __str__(self):
        """

        :return: string representation of the error
        """
        return repr(self.value)


def is_an_oxide_and_no_env_for_O(lse: LightStructureEnvironments) -> bool:
    """

    :param lse: object of the type LightStructureEnvironments
    :return: Boolean if the compound is an oxide and no environments were computed for oxygen
    """
    # TODO: should also check if there are neighbors for every cation!
    for isite, site in enumerate(lse.structure):
        if lse.valences[isite] < 0 and site.species_string != 'O':
            raise ValueError("This is not an oxide. The assessment will be stopped.")
    for isite, site_envs in enumerate(lse.coordination_environments):
        if lse.structure[isite].species_string == 'O':
            if site_envs != None:
                raise ValueError(
                    "Site_envs of anions have been computed. The code has to stop. Use only_cations in compute_structure_environments")
    return True


def get_entropy_from_frequencies(dict_all_environments: dict, max_env=66) -> dict:
    """
    will calulate a relative Shannon entropy from the frequencies of the chemical environments
    :param dict_all_environments: e.g., {'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}
    :param max_env: maximum number of chemical environments, 66
    :return:  a dict for each element that includes the entropy relative to the maximal shannon entropy, e.g. {'Ga':0.5}
    """
    MAX_ENV = max_env
    all = {}
    frequencies = {}
    entropy = {}
    for key, value in dict_all_environments.items():
        all[key] = 0
        frequencies[key] = {}
        entropy[key] = 0.0
        for key1, value1 in Counter(value).items():
            all[key] += value1
        for key1, value1 in Counter(value).items():
            frequencies[key][key1] = float(value1) / float(all[key])
        for key1, value in frequencies[key].items():
            entropy[key] += -math.log(value, 2) * value

    # print(Coordination)
    # something wrong here!
    maxentropy = (-math.log(1.0 / float(MAX_ENV), 2) * 1.0 / float(MAX_ENV)) * MAX_ENV

    perc_entropy = {}
    for key, value in entropy.items():
        perc_entropy[key] = (float(maxentropy) - float(value)) / float(maxentropy)

    return perc_entropy


def get_most_frequent_environment(dict_all_environments: dict) -> dict:
    """
    will calculate the frequency of the most frequent environment
    :param dict_all_environments: e.g., {'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}
    :return: dict that looks the following: {'Ga': 0.5}
    """
    all = {}
    frequencies = {}
    perc_most_frequent = {}
    perc_ready = {}
    for key, value in dict_all_environments.items():
        all[key] = 0
        frequencies[key] = {}

        for key1, value1 in Counter(value).items():
            all[key] += value1
        for key1, value1 in Counter(value).items():
            frequencies[key][key1] = float(value1) / float(all[key])
        perc_most_frequent[key] = OrderedDict(sorted(frequencies[key].items(), key=lambda x: x[1], reverse=True))

        for value1 in perc_most_frequent[key].values():
            perc_ready[key] = value1
            break;
    return perc_ready


def get_mean_CN_from_frequencies(dict_all_environments: dict) -> dict:
    """
    will calculate the mean coordination number based on the dict
    :param dict_all_environments: e.g., {'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}
    :return: dict that looks the following: {'Ga':5}, where 4 is the mean CN of Ga for the given dict
    """
    numbers_ready = {}
    dict_all_CN = {}
    for key, value in dict_all_environments.items():
        if not key in dict_all_CN:
            dict_all_CN[key] = []
        for value2 in value:
            dict_all_CN[key].append(int(value2.split(":")[1]))  # )

    for key, value in dict_all_CN.items():
        numbers_ready[key] = np.mean(value)
    return numbers_ready


class FrequencyEnvironmentPauling1:
    """
    Class to get dicts of the following type: {'Ga': ['O:6', 'O:6', 'T:4', 'T:4']} for each LightStructureEnvironment
    """

    def __init__(self, lse: LightStructureEnvironments):
        """

        :param lse: Objet of the type LightStructureEnvironments
        """
        is_an_oxide_and_no_env_for_O(lse)

        self.output_dict = {}
        for isite, site_envs in enumerate(lse.coordination_environments):
            # identifies cationic sites - only cations have site_envs
            if site_envs != None:
                if len(site_envs) > 0:
                    if not lse.structure[isite].species_string in self.output_dict:
                        self.output_dict[lse.structure[isite].species_string] = []
                    self.output_dict[lse.structure[isite].species_string].append(site_envs[0]['ce_symbol'])

    def get_details(self):
        """
        will give you the results of the computations
        :return: dict similar to {'Ga': ['O:6', 'O:6', 'T:4', 'T:4']}
        """
        return self.output_dict


class Pauling0:
    """
    Class to count the number of cations in a structure
    """

    def __init__(self, lse: LightStructureEnvironments):
        """

        :param lse: LightStructureEnvironments object where oxygen does not have a calculated environment
        """
        is_an_oxide_and_no_env_for_O(lse)
        self.lse = lse

    def get_cations_in_structure(self) -> dict:
        """
        will return a dict with the counted cations
        :return: dict of this type: {'As': 8}
        """
        elements = []
        valences = self.lse.valences
        for isite, site in enumerate(self.lse.structure):
            if valences[isite] >= 0:
                elements.append(self.lse.structure[isite].species_string)
        return Counter(elements)


class Pauling1:
    """
    Class to test Pauling's first rule
    """

    def __init__(self, lse: LightStructureEnvironments, filenameradii="univalent_cat_radii.json", onlylowerlimit=False):
        """
        :param lse: LightStructureEnvironments, only environments for cations should have been calculated
        :param filenameradii: name of the file containing radius ratios
        :param onlylowerlimit: If False, ratio windows are considered
        """
        is_an_oxide_and_no_env_for_O(lse)

        # important parameters
        self.cat_list = {}
        self.cat_valence_list = {}

        self.mat_pauling_fulfilled = 0
        self.env_not = 0
        self.no_env = 0
        self.no_cat = 0

        with open(filenameradii) as dd:
            dict_radii = json.load(dd)
        species = lse.structure.species
        for isite, site_envs in enumerate(lse.coordination_environments):
            # identifies cationic sites - only cations have site_envs
            if site_envs != None:
                if len(site_envs) > 0:
                    try:
                        try:
                            iratio = dict_radii[species[isite].symbol]['ratio_oxide']
                        except:
                            raise ValueError(
                                species[isite].symbol + " not in the Pauling list")
                        if self._first_rule(iratio, site_envs[0], onlylowerlimit=onlylowerlimit):
                            cat_symbol = species[isite].symbol
                            val = lse.valences[isite]

                            if not cat_symbol in self.cat_list:
                                self.cat_list[cat_symbol] = [0, 0]
                            if not cat_symbol in self.cat_valence_list:
                                self.cat_valence_list[cat_symbol] = {}
                            if not val in self.cat_valence_list[cat_symbol]:
                                self.cat_valence_list[cat_symbol][val] = [0, 0]

                            self.cat_list[cat_symbol][0] = \
                                self.cat_list[cat_symbol][0] + 1
                            self.cat_valence_list[cat_symbol][lse.valences[isite]][0] = \
                                self.cat_valence_list[cat_symbol][lse.valences[isite]][0] + 1
                            self.mat_pauling_fulfilled += 1
                        else:
                            cat_symbol = species[isite].symbol
                            val = lse.valences[isite]

                            if not cat_symbol in self.cat_list:
                                self.cat_list[cat_symbol] = [0, 0]
                            if not cat_symbol in self.cat_valence_list:
                                self.cat_valence_list[cat_symbol] = {}
                            if not val in self.cat_valence_list[cat_symbol]:
                                self.cat_valence_list[cat_symbol][val] = [0, 0]

                            self.cat_list[species[isite].symbol][1] = \
                                self.cat_list[species[isite].symbol][1] + 1
                            self.cat_valence_list[species[isite].symbol][lse.valences[isite]][1] = \
                                self.cat_valence_list[species[isite].symbol][lse.valences[isite]][
                                    1] + 1

                            self.env_not += 1
                    except ValueError as err:
                        if err.args[0] == "env not in Pauling book":
                            self.no_env += 1
                        else:
                            self.no_cat += 1

    def get_details(self) -> dict:
        """

        :return: returns a dict with many details on the Pauling rules, looks like this:
        {'cat_dependency': {'As': [4, 4]}, 'cat_val_dependency': {'As': {5: [4, 4]}},
        'Env_fulfilled': 4, 'Env_notfulfilled': 4, 'Env_out_of_list': 0, 'Cat_out_of_list': 0}
        """
        Outputdict = {}
        # first entry in list gives you the environments fulfilling the rule, the second that don't
        Outputdict["cat_dependency"] = self.cat_list
        # first entry in list gives you the environments fulfilling the rule, the second that don't
        Outputdict["cat_val_dependency"] = self.cat_valence_list
        Outputdict["Env_fulfilled"] = self.mat_pauling_fulfilled
        Outputdict["Env_notfulfilled"] = self.env_not
        Outputdict["Env_out_of_list"] = self.no_env
        Outputdict["Cat_out_of_list"] = self.no_cat

        return Outputdict

    def is_fulfilled(self) -> bool:
        """
        tells you if the rule is fulfilled.
        :return: Boolean
        raises TypeError if not all environments are in Pauling's book and if one univalent radii is not present
        """
        if self.no_env != 0 or self.no_cat != 0:
            raise RuleCannotBeAnalyzedError("The first rule cannot be evaluated.")

        if self.env_not == 0 and self.mat_pauling_fulfilled > 0:
            return True
        else:
            return False

    # directly from Pauling's book
    def _predict_env_pauling_window(self, ratio: float) -> list:
        """
        will return a list with the determined environments
        :param ratio: radius ratio
        :return: list with determined environments
        """
        if ratio < 0.225:
            return ['does not exist']
        elif ratio >= 0.225 and ratio < 0.414:
            return ['T:4']  # tetrahedron
        elif ratio >= 0.414 and ratio < 0.592:
            return ['O:6']  # octahedron
        elif ratio >= 0.592 and ratio < 0.645:
            return ['FO:7']  # face-capped octahedron
        elif ratio >= 0.645 and ratio < 0.732:
            return ['SA:8']  # square antiprism
        elif ratio >= 0.732 and ratio < 1.0:
            return ['TT_1:9', 'C:8']  # Tricapped triangular prism (three square-face caps), #cube
        elif ratio >= 1.0:
            return ['C:12']  # Cuboctahedron

    def _predict_env_pauling_lowerlimit(self, ratio: float) -> list:
        """
        will return a list with the determined environments, this time the radius ratio is treated as a lower limit
        :param ratio: radius ratio
        :return: list with determined environments
        """
        ListToReturn = []
        if ratio >= 0.225:
            ListToReturn.append('T:4')  # tetrahedron
        if ratio >= 0.414:
            ListToReturn.append('O:6')  # octahedron
        if ratio >= 0.592:
            ListToReturn.append('FO:7')  # face-capped octahedron
        if ratio >= 0.645:
            ListToReturn.append('SA:8')  # square antiprism
        if ratio >= 0.732:
            ListToReturn.append('TT_1:9')  # Tricapped triangular prism (three square-face caps), #cube
            ListToReturn.append('C:8')
        if ratio >= 1.0:
            ListToReturn.append('C:12')  # Cuboctahedron
        return ListToReturn

    def _first_rule(self, iratio: float, site_env: dict, onlylowerlimit: bool) -> bool:
        """
        will evalute the first Pauling rule
        :param iratio: radius ratio
        :param site_env: dict including the 'ce_symbol' key
        :param onlylowerlimit: will see the radius ratio only as a lower limit
        :return:
        """
        environments = ['T:4', 'O:6', 'FO:7',
                        'SA:8', 'TT_1:9', 'C:8', 'C:12']

        if site_env['ce_symbol'] in environments:
            # only values from Pauling book are considered
            if not onlylowerlimit:
                if str(site_env['ce_symbol']) in self._predict_env_pauling_window(iratio):
                    return True
                else:
                    return False

            else:
                if str(site_env['ce_symbol']) in self._predict_env_pauling_lowerlimit(iratio):
                    return True
                else:
                    return False

        else:
            raise ValueError("env not in Pauling book")


class Pauling2(Pauling0):
    """
        Class to test the electrostatic valence rule
    """

    def __init__(self, lse: LightStructureEnvironments):
        """
        :param lse: LightStructureEnvironment, only cations should have an coordination environment
        """
        is_an_oxide_and_no_env_for_O(lse)

        self.electrostatic_bond_strengths = {}
        for isite, site in enumerate(lse.structure):
            if lse.valences[isite] > 0:
                # print(site)
                # print(lse.coordination_environments[isite])
                nb_set = lse.neighbors_sets[isite][0]
                cn = float(len(nb_set))
                for nb_dict in nb_set.neighb_sites_and_indices:
                    # print(nb_dict)
                    nb_isite = nb_dict['index']
                    if nb_isite not in self.electrostatic_bond_strengths:
                        self.electrostatic_bond_strengths[nb_isite] = []
                    self.electrostatic_bond_strengths[nb_isite].append({'cation_isite': isite,
                                                                        'bond_strength': float(
                                                                            lse.valences[isite]) / cn})
        self.anions_bond_strengths = []
        self.lse = lse
        for isite, site in enumerate(lse.structure):
            if lse.valences[isite] < 0:
                if isite in self.electrostatic_bond_strengths:
                    bond_strengths = [nb['bond_strength']
                                      for nb in self.electrostatic_bond_strengths[isite]]
                    cations_isites = [nb['cation_isite']
                                      for nb in self.electrostatic_bond_strengths[isite]]
                else:
                    bond_strengths = []
                    cations_isites = []
                bond_strengths_sum = sum(bond_strengths)
                self.anions_bond_strengths.append({'anion_isite': isite,
                                                   'bond_strengths': bond_strengths,
                                                   'cations_isites': cations_isites,
                                                   'nominal_oxidation_state': lse.valences[isite],
                                                   'bond_strengths_sum': bond_strengths_sum}
                                                  )

    def is_fulfilled(self, tolerance=1e-2) -> bool:
        """
        Tells you if rule is fulfilled for the whole structure
        :param: tolerance for deviation from 2.
        :return: Boolean, True if the rule is fulfilled
        """
        self.satisfied = True
        ianionsite = 0
        for isite, site in enumerate(self.lse.structure):
            if self.lse.valences[isite] < 0:
                if (abs(self.anions_bond_strengths[ianionsite]['bond_strengths_sum'] - 2.0)) > tolerance:
                    self.satisfied = False

                ianionsite = ianionsite + 1

        return self.satisfied

    def get_details(self, tolerance=10e-2) -> dict:
        """
        Gives you an output dict with information on each anion
        :param tolerance: tolerance for evaluation of fulfillment for each oxygen
        :return: OutputDict with information on each anion
        """
        OutputDict = {}
        OutputDict["bvs_for_each_anion"] = self._get_anions_bvs()
        OutputDict["cations_around_anion"] = self._get_cations_around_anion()
        OutputDict["elementwise_fulfillment"] = self._get_elementwise_fulfillment(tolerance=tolerance)
        OutputDict["cations_in_structure"] = self.get_cations_in_structure()
        return OutputDict

    def _get_anions_bvs(self) -> list:
        """
        get bond valences sums for each anion
        :return: list of bvs
        """
        bvs = []
        ianionsite = 0
        for isite, site in enumerate(self.lse.structure):
            if self.lse.valences[isite] < 0:
                bvs.append(
                    self.anions_bond_strengths[ianionsite]['bond_strengths_sum'])
                ianionsite = ianionsite + 1
        return bvs

    def _get_cations_around_anion(self) -> list:
        """
        :return: list of list with elements around anion
        """
        ianionsite = 0
        elements = []
        for isite, site in enumerate(self.lse.structure):
            if self.lse.valences[isite] < 0:
                elements_around = [self.lse.structure[icat].species_string for icat in
                                   self.anions_bond_strengths[ianionsite]['cations_isites']]
                elements.append(elements_around)
                ianionsite = ianionsite + 1
        return elements

    def _get_elementwise_fulfillment(self, tolerance=10e-2) -> dict:
        """
        will calculate an elementwise fulfillment
        will count cations as fulfilling if they are around an anion that fulfills the rule
        will count cations as not fulfilling if they are around an anion that does not fulfill the rule
        there will be double counting of cations due to that
        :param tolerance:
        :return: dict including information for each cation
        """
        # TODO: think about the tolerance!
        cations_around_anion = self._get_cations_around_anion()
        Elementwise_fulfillment = {}
        for ibvs, bvs in enumerate(self._get_anions_bvs()):
            if abs(bvs - 2.0) <= tolerance:
                for cat in cations_around_anion[ibvs]:
                    if cat not in Elementwise_fulfillment:
                        Elementwise_fulfillment[cat] = [0, 0]
                    Elementwise_fulfillment[cat][0] += 1
            else:
                for cat in cations_around_anion[ibvs]:
                    if cat not in Elementwise_fulfillment:
                        Elementwise_fulfillment[cat] = [0, 0]
                    Elementwise_fulfillment[cat][1] += 1
        return Elementwise_fulfillment


class Pauling2_optimized_environments(Pauling2):
    """
        Class to test the electrostatic valence rule, not only most frequent environment
    """

    def __init__(self, lse: LightStructureEnvironments, perc: float = 0.3):
        """

        :param lse: LightStructureEnvironment, only cations should have an coordination environment
        :param perc: decides which ce are considered. The one with the highest fraction is always considered, then the ones with fraction larger or equal to perc are also considered
        """
        structure = lse.structure
        valences = lse.valences
        neighbors_set = lse.neighbors_sets
        coordination_environments = lse.coordination_environments
        self.lse = lse
        self.structure = structure
        self.valences = valences

        is_an_oxide_and_no_env_for_O(lse)
        equivalent_indices = self._get_symmetry_eq_indices(structure=structure)
        env_indices = self._get_env_indices(coordination_environments=coordination_environments, valences=valences,
                                            perc=perc, equivalent_indices=equivalent_indices)
        combination_array_only_cations = self._get_combinations_cations(equivalent_indices=equivalent_indices,
                                                                        env_indices=env_indices, valences=valences)
        combinations_all_ions_list = self._get_combinations_all_ions(
            combination_array_only_cations=combination_array_only_cations, equivalent_indices=equivalent_indices,
            valences=valences)
        self._calculate_anion_bond_strengths(combinations_all_ions_list=combinations_all_ions_list, structure=structure,
                                             valences=valences, neighbors_sets=neighbors_set)

    def _get_symmetry_eq_indices(self, structure: Structure) -> list:
        """
        will return the indices that are equivalent within a structure
        :param structure: Structure object without symmetry
        :return: list
        """
        spanalyzer = SpacegroupAnalyzer(structure=structure)
        equivalent_indices = spanalyzer.get_symmetrized_structure().equivalent_indices
        return equivalent_indices


    def _get_env_indices(self, coordination_environments: list, valences: list, perc: float,
                         equivalent_indices: list) -> list:
        """
        #gets all indices of relevant environments with a fraction larger perc, if there is none larger perc it will always get at least the environment with the larges fraction
        :param coordination_environments: list of coordination environments from lse
        :param valences: valences
        :param perc: only the ce with the highest fraction are included + the ones that have a fraction equal or larger to perc
        :param equivalent_indices:
        :return:
        """
        env_indices = []
        for number in range(0, len(valences)):
            env_indices.append([])
        for list_index in equivalent_indices:
            if valences[list_index[0]] > 0:
                for same, sameenv in enumerate(list_index):
                    env_indices_local = []
                    env_indices_local.append(0)
                    for ienv, env in enumerate(coordination_environments[list_index[0]][1:], 1):
                        if env['ce_fraction'] >= perc:
                            env_indices_local.append(ienv)
                    env_indices[sameenv] = env_indices_local
            else:
                for same, sameenv in enumerate(list_index):
                    env_indices[sameenv] = [None]
        return env_indices

    def _get_combinations_cations(self, equivalent_indices: list, env_indices: list, valences: list) -> list:
        """
        will get relevant combinations of neighbor_sets for the cations
        :param equivalent_indices: which indices in the structure are equivalent
        :param env_indices: indices of the environments
        :param valences: list of valences
        :return:
        """
        combination_array_only_cations = []
        for list_index in equivalent_indices:
            if valences[list_index[0]] >= 0:
                if len(combination_array_only_cations) == 0:
                    for el in env_indices[list_index[0]]:
                        combination_array_only_cations.append([el])
                else:
                    combination_array_only_cations = self._list_combi_recursive(combination_array_only_cations,
                                                                                env_indices[list_index[0]])
        return combination_array_only_cations

    def _get_combinations_all_ions(self, combination_array_only_cations: list, equivalent_indices: list,
                                   valences: list) -> list:
        """
        will produce a list which includes lists with the numbers of the neighbors_sets for each site
        :param combination_array_only_cations: list which indicates the relevant neighbors_sets for the cations only
        :param equivalent_indices: list of equivalent indices
        :param valences: list of valences
        :return: list with all relevant combinations
        """
        max_sites = len(valences)
        combinations_all_ions_list = []

        for number_combi, combi in enumerate(combination_array_only_cations):
            combinations_all_ions = []
            for isite in range(max_sites):
                combinations_all_ions.append([])
            ilist_index = 0
            for list_index in equivalent_indices:
                if valences[list_index[0]] > 0:
                    for index in list_index:
                        combinations_all_ions[index] = combi[ilist_index]
                    ilist_index += 1
                else:
                    for index in list_index:
                        combinations_all_ions[index] = None

            combinations_all_ions_list.append(combinations_all_ions)
        return combinations_all_ions_list

    def _list_combi_recursive(self, list1: list, list2: list) -> list:
        """
        will combine two lists like this: [[0]] + [0,1] -> [[0,0],[0,1]] and so on
        :param list1: list of list
        :param list2: list to add
        :return: new list
        """
        new_list = []
        for el in list1:
            for el2 in list2:
                new = el.copy()
                new.append(el2)
                new_list.append(new)
        return new_list

    def _calculate_anion_bond_strengths(self, combinations_all_ions_list: list, structure: Structure, valences: list,
                                        neighbors_sets: list):
        """
        will calculate the anion_bond_strengths for all combinations of environments that are given
        :param combinations_all_ions_list: list which includes lists with the numbers of the neighbors_sets for each site
        :param structure: Structure object
        :param valences: list of valences
        :param neighbors_sets: neighbors_sets from the LightStructureEnvironment
        :return:
        """
        self.electrostatic_bond_strengths = []
        self.anions_bond_strengths = []
        for icombination_list, combination_list in enumerate(combinations_all_ions_list):
            self.electrostatic_bond_strengths.append({})
            self.anions_bond_strengths.append([])
            for isite, site in enumerate(structure):
                if valences[isite] >= 0:
                    nb_set = neighbors_sets[isite][combination_list[isite]]
                    cn = float(len(nb_set))
                    for nb_dict in nb_set.neighb_sites_and_indices:
                        nb_isite = nb_dict['index']
                        if nb_isite not in self.electrostatic_bond_strengths[icombination_list]:
                            self.electrostatic_bond_strengths[icombination_list][nb_isite] = []
                        self.electrostatic_bond_strengths[icombination_list][nb_isite].append({'cation_isite': isite,
                                                                                               'bond_strength': float(
                                                                                                   valences[
                                                                                                       isite]) / cn})

            for isite, site in enumerate(structure):
                if valences[isite] < 0:
                    if isite in self.electrostatic_bond_strengths[icombination_list]:
                        bond_strengths = [nb['bond_strength']
                                          for nb in self.electrostatic_bond_strengths[icombination_list][isite]]
                        cations_isites = [nb['cation_isite']
                                          for nb in self.electrostatic_bond_strengths[icombination_list][isite]]
                    else:
                        bond_strengths = []
                        cations_isites = []
                    bond_strengths_sum = sum(bond_strengths)
                    self.anions_bond_strengths[icombination_list].append({'anion_isite': isite,
                                                                          'bond_strengths': bond_strengths,
                                                                          'cations_isites': cations_isites,
                                                                          'nominal_oxidation_state': valences[
                                                                              isite],
                                                                          'bond_strengths_sum': bond_strengths_sum}
                                                                         )

    def is_fulfilled(self, tolerance=1e-2) -> bool:
        """
        Tells you if rule is fulfilled for the whole structure
        :param: tolerance for deviation from 2.
        :return: Boolean, True if the rule is fulfilled
        """
        self.satisfied = [True] * len(self.anions_bond_strengths)
        for ianion, anion_bond_strength in enumerate(self.anions_bond_strengths):
            ianionsite = 0
            for isite, site in enumerate(self.structure):
                if self.valences[isite] < 0:
                    if (abs(anion_bond_strength[ianionsite]['bond_strengths_sum'] - 2.0)) > tolerance:
                        self.satisfied[ianion] = False
                    ianionsite = ianionsite + 1
        return (True in self.satisfied)

    def get_details(self, tolerance=10e-2) -> dict:
        """
        Gives you an output dict with information on each anion
        :param tolerance: tolerance for evaluation of fulfillment for each oxygen
        :return: OutputDict with information on each anion
        """
        # TODO: update, only info of best environments
        # TODO: get best anions around neighbors

        OutputDict = {}
        OutputDict["bvs_for_each_anion"] = self._get_anions_bvs()
        OutputDict["cations_around_anion"] = self._get_cations_around_anion()
        OutputDict["elementwise_fulfillment"] = self._get_elementwise_fulfillment(tolerance=tolerance)
        OutputDict["cations_in_structure"] = self.get_cations_in_structure()
        return OutputDict

    def _get_anions_bvs(self) -> list:
        """
        get bond valences sums for each anion, will return the bvs with the lowest mean absolute deviation from 2.0
        TODO: update
        :return: list of bvs
        """
        # find the bvs for all anion_bond_strenghts and return the one with the lowest mean deviation
        bvs_list = []
        for anion_bond_strength in self.anions_bond_strengths:
            bvs = []
            ianionsite = 0
            for isite, site in enumerate(self.structure):
                if self.valences[isite] < 0:
                    bvs.append(
                        anion_bond_strength[ianionsite]['bond_strengths_sum'])
                    ianionsite = ianionsite + 1
            bvs_list.append(bvs)
        abs_dev = 1000.0
        index = 0
        for ilistel, list_el in enumerate(bvs_list):
            new_abs_dev = np.mean([abs(i - 2.0) for i in list_el])
            if new_abs_dev < abs_dev:
                # keep < because this will prefer more probable environments
                abs_dev = new_abs_dev
                index = ilistel
        return bvs_list[index]

    def _get_index_anion_combi_lowest_bvs(self) -> int:
        """
        will get the index of the anion charges based on the bvs list with the lowest mean absolute deviation
        TODO: update
        :return: index
        """
        # find the bvs for all anion_bond_strengths and return the one with the lowest mean deviation
        bvs_list = []
        for anion_bond_strength in self.anions_bond_strengths:
            bvs = []
            ianionsite = 0
            for isite, site in enumerate(self.lse.structure):
                if self.lse.valences[isite] < 0:
                    bvs.append(
                        anion_bond_strength[ianionsite]['bond_strengths_sum'])
                    ianionsite = ianionsite + 1
            bvs_list.append(bvs)

        abs_dev = 1000.0
        index = 0
        for ilistel, list_el in enumerate(bvs_list):
            new_abs_dev = np.mean([abs(i - 2.0) for i in list_el])
            # print(new_abs_dev)
            if new_abs_dev < abs_dev:
                # keep < because this will prefer more probable environments
                abs_dev = new_abs_dev
                # print(abs_dev)
                index = ilistel
        return index

    def _get_cations_around_anion(self) -> list:
        """
        :return: list of list with elements around anion
        """

        # get correct anions around
        ianionsite = 0
        elements = []
        for isite, site in enumerate(self.structure):
            if self.valences[isite] < 0:
                elements_around = [self.structure[icat].species_string for icat in
                                   self.anions_bond_strengths[self._get_index_anion_combi_lowest_bvs()][ianionsite][
                                       'cations_isites']]
                elements.append(elements_around)
                ianionsite = ianionsite + 1
        return elements

    def _get_elementwise_fulfillment(self, tolerance: float = 10e-2) -> dict:
        """
        will calculate an elementwise fulfillment
        will count cations as fulfilling if they are around an anion that fulfills the rule
        will count cations as not fulfilling if they are around an anion that does not fulfill the rule
        there will be double counting of cations due to that
        :param tolerance:
        :return: dict including information for each cation
        """
        cations_around_anion = self._get_cations_around_anion()
        Elementwise_fulfillment = {}
        for ibvs, bvs in enumerate(self._get_anions_bvs()):
            if abs(bvs - 2.0) <= tolerance:
                for cat in cations_around_anion[ibvs]:
                    if cat not in Elementwise_fulfillment:
                        Elementwise_fulfillment[cat] = [0, 0]
                    Elementwise_fulfillment[cat][0] += 1
            else:
                for cat in cations_around_anion[ibvs]:
                    if cat not in Elementwise_fulfillment:
                        Elementwise_fulfillment[cat] = [0, 0]
                    Elementwise_fulfillment[cat][1] += 1
        return Elementwise_fulfillment


class PaulingConnection:
    """Class that can analyze connections of polyhedra, will be used for rules 3 to 5"""

    def __init__(self, DISTANCE: float):
        """

        :param DISTANCE: maximum distance between the cations that is considered
        """
        self.DISTANCE = DISTANCE

    def _is_cationic_site(self, isite: int, valences: list) -> bool:
        """
        tells you if element on site is a cation
        :param isite: number of the site that is relevant
        :param valences: list of valences in order of sites
        :return: Boolean if the element on the site is a cation
        """
        if valences[isite] >= 0:
            return True
        else:
            return False

    def _get_oxygen_neighbors(self, lse: LightStructureEnvironments, site: PeriodicSite, r: float, CN: int,
                              ceindex: int):
        """
        will return a list of oxygen neighbor sites
        :param lse: LightStructureEnvironment
        :param site: relevant site
        :param r: radius that should be considered for the analysis of the neighbors
        :param CN: coordination number of this cationic site
        :param ceindex: index of the coordination number at hand (for the most important one, usually 0)
        :return: list of oxygen neighbor sites
        """
        struct = lse.structure
        all_neighbors = struct.get_neighbors(site=site, r=r)
        sites_and_indices = lse.neighbors_sets[self._get_site_index(
            site, struct)][ceindex].neighb_sites_and_indices
        indices = [x['index'] for x in sites_and_indices]
        sorted_all = sorted(all_neighbors, key=lambda x: x[1])
        counter = 0
        oxygen_neighbors = []
        for i in sorted_all:
            if self._get_site_index(i[0], struct) in indices and counter < CN:
                oxygen_neighbors.append(i[0])
                counter = counter + 1
                indices.remove(self._get_site_index(i[0], struct))
            if counter == CN:
                break
        return oxygen_neighbors

    def _get_cation_neighbors(self, struct: Structure, site: PeriodicSite, r: float, valences: list) -> list:
        """
        will return a list of neighboring cationic neighbor sites
        get cationic neighbors of a cation
        :param struct: Structure Object
        :param site: relevant cationic site
        :param r: radius
        :param valences: list of valences
        :return: list of neighboring sites
        """
        all_neighbors = struct.get_neighbors(
            site=site, r=r, include_index=False)
        cation_neighbors = []
        for i in all_neighbors:
            if self._is_cationic_site(self._get_site_index(i[0], struct), valences):
                cation_neighbors.append(i[0])
        return cation_neighbors

    def _get_site_index(self, insite: PeriodicSite, struct: Structure) -> int:
        """
        will return index of a site
        :param insite: PeriodicSite
        :param struct: Structure
        :return: integer indicating the site index
        """
        sites = struct.sites
        for isite, site in enumerate(sites):
            if insite.is_periodic_image(site):
                return isite


class Pauling3and4(PaulingConnection):
    """
    Class to analyse third and fourth rule
    """

    def __init__(self):
        """
        will overwrite the PaulingConnection init
        """
        pass

    def newsetup(self, lse: LightStructureEnvironments, filename=None, save_to_file=True,
                 foldername='ThirdRuleAnalysisConnections', distance=8.0):
        """
            will setup everything from scratch, should only be used in the beginning of an analysis ->slow
            :param lse: LightStructureEnvironment
            :param save_to_file: Boolean
            :param filename: usually string, for example "file1.json"
            :param foldername: name of the folder you would like to save the file in
            :param distance: float giving the distances of cations that is considered
        """
        self.lse = lse
        is_an_oxide_and_no_env_for_O(lse)
        super().__init__(DISTANCE=distance)
        if filename is None:
            save_to_file == False
        if save_to_file:
            if not os.path.isdir(foldername) and (foldername != ''):
                os.mkdir(foldername)
        DISTANCE = self.DISTANCE

        # get all relevant sites
        allsites = self._get_allsites(lse, DISTANCE)
        # search all pairs of sites that fulfill the DISTANCE criterion
        pairedsites = self._get_pairedsites(allsites, DISTANCE)
        self.PolyhedronDict = self._get_connections(pairedsites, lse, DISTANCE)
        if save_to_file:
            with open(os.path.join(foldername, filename), 'w') as file:
                json.dump(self.PolyhedronDict, file)

    def from_file(self, filename: str, foldername='ThirdRuleAnalysisConnections'):
        """
        will setup everything from a json file
        :param filename: "example.json"
        :param foldername: folder, in which the file is located
        :return:
        """
        with open(os.path.join(foldername, filename)) as file:
            self.PolyhedronDict = json.load(file)

    def _get_connections_Pauling3_and_4(self, maxCN=None) -> dict:
        """
        internal method to get the numbers of connections
        :param maxCN: maximum Coordination number that is considered
        :return: returns connections as a dict
        """
        inputdict = self.PolyhedronDict
        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        not_connected = 0
        corner = 0
        edge = 0
        face = 0
        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if maxCN is None:
                if herepolyhedra[numberpolyhedra] == 0:
                    not_connected += 1
                if herepolyhedra[numberpolyhedra] == 1:
                    corner += 1
                if herepolyhedra[numberpolyhedra] == 2:
                    edge += 1
                if herepolyhedra[numberpolyhedra] >= 3:
                    face += 1
            else:
                if (iinfo['CN'][0] <= maxCN) and (iinfo['CN'][1] <= maxCN):
                    if herepolyhedra[numberpolyhedra] == 0:
                        not_connected += 1
                    if herepolyhedra[numberpolyhedra] == 1:
                        corner += 1
                    if herepolyhedra[numberpolyhedra] == 2:
                        edge += 1
                    if herepolyhedra[numberpolyhedra] >= 3:
                        face += 1

            numberpolyhedra += 1

        return {"no": not_connected, "corner": corner, "edge": edge, "face": face}

    def _get_number_connected_oxygens_and_CN(self, index: int, index2: int, sitetupel: tuple,
                                             lse: LightStructureEnvironments) -> tuple:
        """
        tells you how many oxygen atoms are shared by two polyhedra
        :param index: site index of first cation site
        :param index2: site index of second cation site
        :param sitetupel: tuple of sites that is investigated
        :param lse: LightStructureEnvironment
        :return: returns number of connected oxygens, coordination number of first cation and second cation
        """
        maxce = max(
            lse.coordination_environments[index], key=lambda x: x['ce_fraction'])
        numberce = lse.coordination_environments[index].index(maxce)
        CN = int(maxce['ce_symbol'].split(":")[1])
        oxygensites = self._get_oxygen_neighbors(
            lse, sitetupel[0], 8, CN, numberce)

        maxce2 = max(
            lse.coordination_environments[index2], key=lambda x: x['ce_fraction'])
        numberce2 = lse.coordination_environments[index].index(maxce)
        CN2 = int(maxce2['ce_symbol'].split(":")[1])

        oxygensites2 = self._get_oxygen_neighbors(
            lse, sitetupel[1], 8, CN2, numberce2)
        intersectOx = list(set(oxygensites).intersection(oxygensites2))
        numberconnect = len(intersectOx)
        return numberconnect, CN, CN2

    def _get_cationvalences(self, valences: list) -> list:
        """
        get valences of cation
        :param valences: valences of all sites
        :return: returns list of valences of cationic sites
        """
        cationvalences = []
        for i, val in enumerate(valences):
            if self._is_cationic_site(i, valences):
                cationvalences.append(val)
        return cationvalences

    def _get_cationCN(self, lse: LightStructureEnvironments):
        """
        :param lse: LightStructureEnvironment
        :return: returns coordination numbers of cations as a list
        """
        sites = lse.structure.sites
        CNlist = []
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=lse.valences):
                maxce = max(
                    lse.coordination_environments[isite], key=lambda x: x['ce_fraction'])
                CN = int(maxce['ce_symbol'].split(":")[1])
                CNlist.append(CN)
        return CNlist

    def _get_connections(self, pairedsites: list, lse: LightStructureEnvironments, DISTANCE: float) -> dict:
        """
        will give an output dict with information on the connection of polyhedra
        :param pairedsites: paired sites
        :param lse: LightStructureEnvironments
        :param DISTANCE: distance for the connections
        :return: Outputdict with many information
        """
        struct = lse.structure

        valences = lse.valences

        # get list of cation valences
        cationvalences = self._get_cationvalences(valences)
        # is the valence equal
        samevalences = (cationvalences.count(
            cationvalences[0]) == len(cationvalences))
        # get list of coordination numbers for each cation
        CNlist = self._get_cationCN(lse)
        sameCN = (CNlist.count(CNlist[0]) == len(CNlist))

        # do something else
        additionalinfo = []
        herepolyhedra = []
        for sitetupel in pairedsites:
            mydict = {}
            mydict['valences'] = []
            mydict['CN'] = []
            mydict['cations'] = []
            index = self._get_site_index(sitetupel[0], struct)
            index2 = self._get_site_index(sitetupel[1], struct)
            numberconnect, CN, CN2 = self._get_number_connected_oxygens_and_CN(
                index, index2, sitetupel, lse)

            mydict['valences'] = [valences[index], valences[index2]]

            mydict['CN'] = [CN, CN2]

            mydict['distance'] = np.linalg.norm(
                sitetupel[0].coords - sitetupel[1].coords)
            cat1 = sitetupel[0].species_string
            cat2 = sitetupel[1].species_string
            mydict['cations'] = [cat1, cat2]
            additionalinfo.append(mydict)
            herepolyhedra.append(numberconnect)

        numberofPolyhedra = len(herepolyhedra)
        nofPoly_notconnected = len([x for x in herepolyhedra if x == 0])
        nofPoly_corner = len([x for x in herepolyhedra if x == 1])
        noofPoly_edge = len([x for x in herepolyhedra if x == 2])
        noofPoly_face = len([x for x in herepolyhedra if x >= 3])

        outputdict = {}
        outputdict['MaxDistance'] = DISTANCE
        outputdict['Polyhedra'] = numberofPolyhedra
        outputdict['Not'] = nofPoly_notconnected
        outputdict['Corner'] = nofPoly_corner
        outputdict['Edge'] = noofPoly_edge
        outputdict['Face'] = noofPoly_face
        outputdict['PolyConnect'] = herepolyhedra
        outputdict['Additional'] = additionalinfo
        outputdict['cationvalences'] = cationvalences
        outputdict['samevalences'] = samevalences
        outputdict['CNlist'] = CNlist
        outputdict['sameCN'] = sameCN

        return outputdict

    def _get_pairedsites(self, allsites: list, DISTANCE: float) -> list:
        """
        will pair sites that have a certain distance
        :param allsites: all sites that will be considered
        :param DISTANCE: only sites having a distance smaller or equal to this distance will be considered
        :return: list of paired sites in list
        """
        startendpoints = []
        pairedsites = []
        for i in range(0, len(allsites)):
            for j in range(i, len(allsites)):
                #
                fraccoorda = allsites[i].frac_coords.copy()
                fraccoordb = allsites[j].frac_coords.copy()

                if not (np.linalg.norm(allsites[i].coords - allsites[j].coords) > DISTANCE or (
                        fraccoorda[0] == fraccoordb[0] and fraccoorda[1] == fraccoordb[1] and fraccoorda[2] ==
                        fraccoordb[2])):

                    # make sure that at least one coordinate lies in the first cell [0:1)
                    while fraccoorda[0] < 0.0 or fraccoorda[1] < 0.0 or fraccoorda[2] < 0.0 or fraccoordb[
                        0] < 0.0 or fraccoordb[1] < 0.0 or fraccoordb[2] < 0.0:
                        # this code part might be redundant
                        for ij in range(0, 3):
                            if fraccoorda[ij] < 0.0 or fraccoordb[ij] < 0.0:
                                fraccoorda[ij] = fraccoorda[ij] + 1.0
                                fraccoordb[ij] = fraccoordb[ij] + 1.0
                    while not ((fraccoorda[0] >= 0.0 and fraccoordb[0] >= 0.0)
                               and (fraccoorda[0] < 1.0 or fraccoordb[0] < 1.0)
                               and (fraccoorda[1] >= 0.0 and fraccoordb[1] >= 0.0)
                               and (fraccoorda[1] < 1.0 or fraccoordb[1] < 1.0)
                               and (fraccoorda[2] >= 0.0 and fraccoordb[2] >= 0.0)
                               and (fraccoorda[2] < 1.0 or fraccoordb[2] < 1.0)):
                        for number in range(0, 3):
                            if ((fraccoorda[number] >= 1.0) and fraccoordb[number] >= 1.0):
                                # print('test2')
                                fraccoorda[number] = fraccoorda[number] - 1.0
                                fraccoordb[number] = fraccoordb[number] - 1.0

                    # Sort to decide whether we already considered this pair
                    sort = sorted([fraccoorda, fraccoordb],
                                  key=lambda x: (x[0], x[1], x[2]))
                    contained = False
                    for istartendpoints in startendpoints:
                        if np.allclose(sort, istartendpoints):
                            contained = True
                            break
                    if not contained:
                        startendpoints.append(sort)
                        pairedsites.append([allsites[i], allsites[j]])
        return pairedsites

    def _get_allsites(self, lse: LightStructureEnvironments, DISTANCE: float) -> list:
        """
        will help investigate the connections
        :param lse: LightStructureEnvironment
        :param DISTANCE: distance between cations
        :return: returns relevant sites to investigate the connections as a list
        """

        sites = lse.structure.sites
        valences = lse.valences
        struct = lse.structure
        allsites = []
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=valences):
                # get infos from site
                atoms_n_occu = site.species_string
                lattice = site.lattice
                coords = site.frac_coords

                # get number of cells that have to be considered to find neighbors with a CERTAIN distance
                # inspired by Lattice class in pymatgen.core.lattice
                recp_len = np.array(struct.lattice.reciprocal_lattice.abc) / (2 * np.pi)
                nmax = float(DISTANCE) * recp_len + 0.01

                # +1 due to the fact that also neighbors of an atom on (1,1,1) have to be considered
                ceila = int(np.ceil(nmax[0]) + 1)
                ceilb = int(np.ceil(nmax[1]) + 1)
                ceilc = int(np.ceil(nmax[2]) + 1)

                for i in range(0, ceila):
                    for j in range(0, ceilb):
                        for k in range(0, ceilc):
                            addcoord = [float(i), float(j), float(k)]
                            newcoords = coords + addcoord
                            allsites.append(PeriodicSite(
                                atoms_n_occu, newcoords, lattice, to_unit_cell=False))

        return allsites

    def get_connections(self, maximumCN=None) -> dict:
        """
        Gives connections
        :return: Outputdict with number of connections
        """
        OutputDict = self._get_connections_Pauling3_and_4(maxCN=maximumCN)
        return OutputDict


class Pauling3(Pauling3and4):
    """
    Class to evaluate the third rule
    """

    def is_fulfilled(self, maximumCN=None) -> bool:
        """
        tells you if third rule is fulfilled (i.e. no connections via faces!)
        :param maximumCN: gives the maximal CN of cations that are considered
        :return: True if the rule is fulfilled (i.e. no connections via faces)
        """
        connections = self._get_connections_Pauling3_and_4(maxCN=maximumCN)

        if connections["face"] == 0:
            return True
        else:
            return False

    def get_details(self, maximumCN=None) -> dict:
        """
        returns a dict that looks the following {"Element string": "Valence as an int": {"no": number1, "corner": number2, "edge": number3, "face" number4} indicating the connections.}
        :param maximumCN: maximum CN that is considered
        :return: Outputdict with information on the types of connetions of pairs of polyhedra and with information on the valences, species that are connected via corners, edges, and faces or not connected

        """
        Outputdict = self._postevaluation3rdrule_evaluate_elementdependency(maxCN=maximumCN)
        seconddict = self._get_connections_Pauling3_and_4(maxCN=maximumCN)
        seconddict["species"] = Outputdict

        return seconddict

    def _postevaluation3rdrule_evaluate_elementdependency(self, maxCN=None) -> dict:
        """
        evaluates element and valence dependency of connected pairs of polyhedra
        :param maxCN:
        :return: a Outputdict that will be evaluated by other methods
        """
        inputdict = self.PolyhedronDict
        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        # Dict that includes cations sharing no connections, corners, edges, and faces
        # This can be used to evaluate the results per polyhedron pair
        Dict_ElementDependency = {}
        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if maxCN is None:
                if iinfo['cations'][0] not in Dict_ElementDependency:
                    Dict_ElementDependency[iinfo['cations'][0]] = {}
                if iinfo['cations'][1] not in Dict_ElementDependency:
                    Dict_ElementDependency[iinfo['cations'][1]] = {}

                if iinfo['valences'][0] not in Dict_ElementDependency[iinfo['cations'][0]]:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]] = {'no': 0, 'corner': 0,
                                                                                         'edge': 0, 'face': 0}

                if iinfo['valences'][1] not in Dict_ElementDependency[iinfo['cations'][1]]:
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]] = {'no': 0, 'corner': 0,
                                                                                         'edge': 0, 'face': 0}

                if herepolyhedra[numberpolyhedra] == 0:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['no'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['no'] += 1

                elif herepolyhedra[numberpolyhedra] == 1:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['corner'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['corner'] += 1


                elif herepolyhedra[numberpolyhedra] == 2:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['edge'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['edge'] += 1


                elif herepolyhedra[numberpolyhedra] > 2:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['face'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['face'] += 1

            else:
                if iinfo["CN"][0] <= maxCN and iinfo["CN"][1] <= maxCN:
                    if iinfo['cations'][0] not in Dict_ElementDependency:
                        Dict_ElementDependency[iinfo['cations'][0]] = {}
                    if iinfo['cations'][1] not in Dict_ElementDependency:
                        Dict_ElementDependency[iinfo['cations'][1]] = {}

                    if iinfo['valences'][0] not in Dict_ElementDependency[iinfo['cations'][0]]:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]] = {'no': 0, 'corner': 0,
                                                                                             'edge': 0, 'face': 0}

                    if iinfo['valences'][1] not in Dict_ElementDependency[iinfo['cations'][1]]:
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]] = {'no': 0, 'corner': 0,
                                                                                             'edge': 0, 'face': 0}

                    if herepolyhedra[numberpolyhedra] == 0:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['no'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['no'] += 1

                    elif herepolyhedra[numberpolyhedra] == 1:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['corner'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['corner'] += 1


                    elif herepolyhedra[numberpolyhedra] == 2:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['edge'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['edge'] += 1


                    elif herepolyhedra[numberpolyhedra] > 2:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['face'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['face'] += 1

            numberpolyhedra = numberpolyhedra + 1

        return Dict_ElementDependency


class Pauling4(Pauling3and4):
    """
    Class to evaluate the Fourth Rule
    """

    def is_fulfilled(self) -> bool:
        """
        tells you if polyhedra of cations with highest valence and smallest CN don't show any connections within the structure
        structure has to have cations with different valences and coordination numbers
        :return: Boolean
        :raises: RuleCannotBeAnalyzedError if structure cannot will be accessed for this test (no cation that differ in valence and CN)
        """
        if not self._is_candidate_4thrule():
            raise RuleCannotBeAnalyzedError("The fourth rule cannot be evaluated")

        maxval = max(self.PolyhedronDict['cationvalences'])
        minCN = min(self.PolyhedronDict['CNlist'])
        outputdict = self._postevaluation4thruleperpolyhedron_only_withoutproduct(minCN, minCN, maxval, maxval)
        if outputdict['corner'] != 0 or outputdict['edge'] != 0 or outputdict['face'] != 0:
            return False
        else:
            return True

    def get_details(self) -> dict:
        """
        gives you number of connections as a function of the coordination numbers and valences of the polyhedra pairs
        :return: Dict of the following form: {"val1:valence1": {"val2:valence2": {"CN1:CN1": {"CN2:CN2": {"no": number1, "corner": number2, "edge": number3, "face" number4}}}}} indicating the connections depending on valences and CN
        """
        if not self._is_candidate_4thrule():
            raise RuleCannotBeAnalyzedError()
        return self._postevaluation4thrule()

    def _is_candidate_4thrule(self) -> bool:
        inputdict = self.PolyhedronDict
        samevalences = inputdict['samevalences']
        sameCN = inputdict['sameCN']
        if (not samevalences) or (not sameCN):
            return True
        else:
            return False

    # Add additional evaluation of 4th rule that does not depend on valence or CN product but on 2 valences, 2 CN
    def _postevaluation4thrule(self) -> dict:
        """
        :return: more information connected pairs of polyhedra
        """
        # TODO: make additional tests here
        # TODO: there might be a small bug here
        inputdict = self.PolyhedronDict

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        Outputdict = {}
        Elementwise = {}
        numberpolyhedra = 0
        for iinfo in additionalinfo:

            if not iinfo['cations'][0] in Elementwise:
                Elementwise[iinfo['cations'][0]] = {}
            if not iinfo['cations'][1] in Elementwise:
                Elementwise[iinfo['cations'][1]] = {}

            # here: consider elements:
            # symmetry of valence and CN have to be considered
            # have to think about symmetry again!!!!
            if not ("val1:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][0]]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])] = {}
            if not ("val1:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][0]]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][0]][
                "val1:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][0]][
                "val1:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])] = {}

            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])] = {}

            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}

            if herepolyhedra[numberpolyhedra] == 0:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
            elif herepolyhedra[numberpolyhedra] == 1:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1

            elif herepolyhedra[numberpolyhedra] == 2:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1

            elif herepolyhedra[numberpolyhedra] > 2:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1

            # now, let's consider the second cation:
            if not ("val1:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][1]]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])] = {}
            if not ("val1:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][1]]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][1]][
                "val1:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][1]][
                "val1:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])] = {}

            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])] = {}

            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}

            if herepolyhedra[numberpolyhedra] == 0:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
            elif herepolyhedra[numberpolyhedra] == 1:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1

            elif herepolyhedra[numberpolyhedra] == 2:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1

            elif herepolyhedra[numberpolyhedra] > 2:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1

            # symmetry of valence and CN have to be considered
            if not ("val1:" + str(iinfo['valences'][0])) in Outputdict:
                Outputdict["val1:" + str(iinfo['valences'][0])] = {}
            if not ("val1:" + str(iinfo['valences'][1])) in Outputdict:
                Outputdict["val1:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][1])) in Outputdict["val1:" + str(iinfo['valences'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][0])) in Outputdict["val1:" + str(iinfo['valences'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])] = {}

            if not ("CN1:" + str(iinfo['CN'][0])) in Outputdict["val1:" + str(iinfo['valences'][0])][
                "val2:" + str(iinfo['valences'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][0])) in Outputdict["val1:" + str(iinfo['valences'][1])][
                "val2:" + str(iinfo['valences'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in Outputdict["val1:" + str(iinfo['valences'][0])][
                "val2:" + str(iinfo['valences'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in Outputdict["val1:" + str(iinfo['valences'][1])][
                "val2:" + str(iinfo['valences'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])] = {}

            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}

            if herepolyhedra[numberpolyhedra] == 0:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
            elif herepolyhedra[numberpolyhedra] == 1:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1

            elif herepolyhedra[numberpolyhedra] == 2:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1

            elif herepolyhedra[numberpolyhedra] > 2:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1

            numberpolyhedra = numberpolyhedra + 1

        Outputdict['maxval'] = max(self.PolyhedronDict['cationvalences'])
        Outputdict['minCN'] = min(self.PolyhedronDict['CNlist'])

        Outputdict['elementwise'] = Elementwise

        return Outputdict

    def _postevaluation4thruleperpolyhedron_only_withoutproduct(self, CN1: int, CN2: int, val1: int, val2: int) -> dict:
        """
        That is the one I used for the evaluation
        :param CN1: CN1 of an atom in pair of polyhedra
        :param CN2: CN2 of the other atom in pair of polyhedra
        :param val1: valence 1 of an atom in pair of polyhedra
        :param val2: valence 2 of an atom in pair of polyhedra
        :return: more information connected pairs of polyhedra
        """
        inputdict = self.PolyhedronDict

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        notconnected = 0
        corner = 0
        edge = 0
        face = 0

        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if ((iinfo['valences'][0] == val1) and (iinfo['valences'][1] == val2)) or (
                    (iinfo['valences'][1] == val1) and (iinfo['valences'][0] == val2)):
                if ((iinfo['CN'][0] == CN1 and iinfo['CN'][1] == CN2) or (
                        iinfo['CN'][0] == CN2 and iinfo['CN'][1] == CN1)):
                    if (iinfo['valences'][0] == val1 and iinfo['CN'][0] == CN1) or (
                            iinfo['valences'][1] == val1 and iinfo['CN'][1] == CN1) or (
                            iinfo['valences'][0] == val2 and iinfo['CN'][0] == CN2) or (
                            iinfo['valences'][1] == val2 and iinfo['CN'][1] == CN2):
                        if herepolyhedra[numberpolyhedra] == 0:
                            notconnected = notconnected + 1
                        elif herepolyhedra[numberpolyhedra] == 1:
                            corner = corner + 1
                        elif herepolyhedra[numberpolyhedra] == 2:
                            edge = edge + 1
                        elif herepolyhedra[numberpolyhedra] > 2:
                            face = face + 1

            numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['no'] = notconnected
        Thirdruledict['corner'] = corner
        Thirdruledict['edge'] = edge
        Thirdruledict['face'] = face

        return Thirdruledict


class Pauling5(PaulingConnection):
    """
    Class to evaluate the fifth rule
    """

    def __init__(self):
        pass

    def newsetup(self, lse: LightStructureEnvironments, filename=None, save_to_file=True,
                 foldername="FifthsRuleAnalysis", distance=8.0):
        """
        collects the coordination numbers, coordination environments and the number of connections via corners, edges and faces of each of the polyhedra for each of the cations

        :param lse: LightStructureEnvironment
        :param filename: "test.json", name under which the files are saved
        :param save_to_file: if this is True, the files will be saved
        :param foldername: name of the folder that the files will be saved in
        :param distance: maximum distance that is considered in evaluation of connections
        :return:
        """

        is_an_oxide_and_no_env_for_O(lse)
        super().__init__(DISTANCE=distance)

        if filename is None:
            save_to_file == False
        if save_to_file:
            if not os.path.isdir(foldername) and (foldername != ''):
                os.mkdir(foldername)
        struct = lse.structure
        sites = struct.sites
        catid = []
        catenv = []
        connection_corners = []
        connection_edges = []
        connection_faces = []
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=lse.valences):
                corner = 0
                edge = 0
                face = 0
                catid.append([site.species_string, lse.valences[isite]])
                maxce = max(
                    lse.coordination_environments[isite], key=lambda x: x['ce_fraction'])
                numberce = lse.coordination_environments[isite].index(maxce)
                CN = int(maxce['ce_symbol'].split(":")[1])
                oxygensites = self._get_oxygen_neighbors(
                    lse, site, 8, CN, numberce)
                catenv.append(maxce['ce_symbol'])
                cation_neighbors = self._get_cation_neighbors(
                    struct, site, r=self.DISTANCE, valences=lse.valences)
                for neigh in cation_neighbors:
                    indexneigh = self._get_site_index(neigh, struct)
                    maxce2 = max(
                        lse.coordination_environments[indexneigh], key=lambda x: x['ce_fraction'])
                    numberce2 = lse.coordination_environments[indexneigh].index(
                        maxce2)
                    CN2 = int(maxce2['ce_symbol'].split(":")[1])
                    oxygensites2 = self._get_oxygen_neighbors(lse, neigh, 8, CN2, numberce2)
                    intersectOx = list(
                        set(oxygensites).intersection(oxygensites2))
                    numberconnect = len(intersectOx)
                    if numberconnect == 1:
                        corner = corner + 1
                    if numberconnect == 2:
                        edge = edge + 1
                    if numberconnect >= 3:
                        face = face + 1
                connection_corners.append(corner)
                connection_edges.append(edge)
                connection_faces.append(face)

        # unique cations with valences
        uniquecat = []
        for cat in catid:
            if cat not in uniquecat:
                uniquecat.append(cat)

        outputdict = {}
        outputdict['connection_corners'] = connection_corners
        outputdict['connection_edges'] = connection_edges
        outputdict['connection_faces'] = connection_faces
        outputdict['catenv'] = catenv
        outputdict['catid'] = catid
        outputdict['uniquecat'] = uniquecat

        if save_to_file:
            with open(os.path.join(foldername, filename), 'w') as file:
                json.dump(outputdict, file)

        self.FifthRuleDict = outputdict

    def from_file(self, filename: str, foldername: str):
        """
        setup of the connections from a saved file
        :param filename: file name under which the file is saved
        :param foldername: folder, in which the files are
        :return:
        """
        with open(os.path.join(foldername, filename)) as file:
            self.FifthRuleDict = json.load(file)

    def is_fulfilled(self, options="CN", leave_out_list=[]) -> bool:
        """
        tests if all chemically equivalent cations (same element, same valence) have the same CN, or environment, or environment and number of connections
        if there is only one chemically equivalent cation, the method raises a RuleCannotBeAnalyzedError
        if the wrong option is used, a ValueError is raised
        :param options: can be "CN", "env", or "env+nconnections"
        :param leave_out_list: list of cations (str) that will be left out from analysis
        :return: Boolean
        """
        # test if rule can be evaluated (more than one chemically equivalent cation!)
        if not self._is_candidate_5thrule(leave_out_list=leave_out_list):
            raise RuleCannotBeAnalyzedError("5th Rule cannot be evaluated")

        details = self.get_details(options=options, leave_out_list=leave_out_list)

        for key, items in details.items():
            if items['not_fulfilled'] > 0:
                return False
        return True

    def get_details(self, options='CN', leave_out_list=[]) -> dict:
        """
        will return a dict with information to assess this rule
        :param options:  can be "CN", "env", or "env+nconnections"
        :param leave_out_list: list of cations (str) that will be left out from analysis
        :return: dict with information on the fulfillment of the rule
        """
        if not self._is_candidate_5thrule(leave_out_list=leave_out_list):
            raise RuleCannotBeAnalyzedError("5th Rule cannot be evaluated")

        outputdict = self._postevaluation5thrule_elementdependency()
        output = {}

        if options == 'CN':
            for cat in outputdict['exceptionsCN']:
                if not cat[0] in leave_out_list:
                    if not cat[0] in output:
                        output[cat[0]] = {}
                        output[cat[0]]['not_fulfilled'] = 0
                        output[cat[0]]['fulfilled'] = 0
                    output[cat[0]]['not_fulfilled'] += 1
            for cat in outputdict['fulfillingCN']:
                if not cat[0] in leave_out_list:
                    if not cat[0] in output:
                        output[cat[0]] = {}
                        output[cat[0]]['not_fulfilled'] = 0
                        output[cat[0]]['fulfilled'] = 0
                    output[cat[0]]['fulfilled'] += 1
            return output
        elif options == 'env':
            for cat in outputdict['exceptionsenvs']:
                if not cat[0] in leave_out_list:
                    if not cat[0] in output:
                        output[cat[0]] = {}
                        output[cat[0]]['not_fulfilled'] = 0
                        output[cat[0]]['fulfilled'] = 0
                    output[cat[0]]['not_fulfilled'] += 1
            for cat in outputdict['fulfillingenvs']:
                if not cat[0] in leave_out_list:
                    if not cat[0] in output:
                        output[cat[0]] = {}
                        output[cat[0]]['not_fulfilled'] = 0
                        output[cat[0]]['fulfilled'] = 0
                    output[cat[0]]['fulfilled'] += 1
            return output
        elif options == 'env+nconnections':
            for cat in outputdict['exceptionsenvs_connections']:
                if not cat[0] in leave_out_list:
                    if not cat[0] in output:
                        output[cat[0]] = {}
                        output[cat[0]]['not_fulfilled'] = 0
                        output[cat[0]]['fulfilled'] = 0
                    output[cat[0]]['not_fulfilled'] += 1
            for cat in outputdict['fulfillingenvs_connections']:
                if not cat[0] in leave_out_list:
                    if not cat[0] in output:
                        output[cat[0]] = {}
                        output[cat[0]]['not_fulfilled'] = 0
                        output[cat[0]]['fulfilled'] = 0
                    output[cat[0]]['fulfilled'] += 1
            return output
        else:
            raise ValueError("Wrong option")

    def _is_candidate_5thrule(self, leave_out_list=[]) -> bool:
        """
        tells you if 5th rule can be evaluated
        :param leave_out_list: will remove those structures that only include cations from the elements in list
        :return: Boolean
        """

        catid = self.FifthRuleDict['catid']
        uniquecat = self.FifthRuleDict['uniquecat']
        if leave_out_list == []:
            if (len(uniquecat) != len(catid)):
                return True
            else:
                return False
        else:
            newcatid = []
            newuniquecat = []
            for cat in catid:
                if cat[0] not in leave_out_list:
                    newcatid.append(cat)
            for cat in uniquecat:
                if cat[0] not in leave_out_list:
                    newuniquecat.append(cat)
            if (len(newuniquecat) != len(newcatid)):
                return True
            else:
                return False

    def _postevaluation5thrule_elementdependency(self) -> dict:
        """
        gives you information on the element dependency of the rule
        :return: dict with information on the elements
        """

        connection_corners = self.FifthRuleDict['connection_corners']
        connection_edges = self.FifthRuleDict['connection_edges']
        connection_faces = self.FifthRuleDict['connection_faces']
        catenv = self.FifthRuleDict['catenv']
        catid = self.FifthRuleDict['catid']
        uniquecat = self.FifthRuleDict['uniquecat']

        hassamechemenvandsameconnectionnumber = True
        hassameenv = True
        hassameCN = True

        exceptions_CN = []
        exceptions_envs = []
        exceptions_envs_connections = []

        if len(uniquecat) != len(catid):

            for icat, cat in enumerate(catid):
                for icat2, cat2 in enumerate(catid):
                    if icat2 > icat:
                        if cat == cat2:
                            if not (int(str(catenv[icat].split(":")[1])) == int(
                                    str(catenv[icat2].split(":")[1]))):
                                hassameCN = False
                                if cat not in exceptions_CN:
                                    exceptions_CN.append(cat)
                            if not (str(catenv[icat]) == str(catenv[icat2])):
                                hassameenv = False
                                if cat not in exceptions_envs:
                                    exceptions_envs.append(cat)
                            if not (connection_corners[icat] == connection_corners[icat2] and connection_edges[
                                icat] ==
                                    connection_edges[icat2] and connection_faces[icat] == connection_faces[
                                        icat2] and str(catenv[icat]) == str(catenv[icat2])):
                                # print('here falsch')
                                hassamechemenvandsameconnectionnumber = False

                                if cat not in exceptions_envs_connections:
                                    exceptions_envs_connections.append(cat)

        fulfilling_CN = []
        fulfilling_envs = []
        fulfilling_envs_connections = []
        for unicat in uniquecat:
            if unicat not in exceptions_CN:
                fulfilling_CN.append(unicat)
            if unicat not in exceptions_envs:
                fulfilling_envs.append(unicat)
            if unicat not in exceptions_envs_connections:
                fulfilling_envs_connections.append(unicat)

        Outputdict = {}

        # checks if number of cations is equal to the number of unique cations
        # find a better way for this CNsokay part
        Outputdict['hassameCN'] = hassameCN
        Outputdict['hassameenv'] = hassameenv
        Outputdict['hassameenvadsameconnectionnumber'] = hassamechemenvandsameconnectionnumber
        Outputdict['exceptionsCN'] = exceptions_CN
        Outputdict['exceptionsenvs'] = exceptions_envs
        Outputdict['exceptionsenvs_connections'] = exceptions_envs_connections
        Outputdict['fulfillingCN'] = fulfilling_CN
        Outputdict['fulfillingenvs'] = fulfilling_envs
        Outputdict['fulfillingenvs_connections'] = fulfilling_envs_connections

        return Outputdict

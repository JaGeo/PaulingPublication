from PaulingRules import Pauling0, Pauling1, Pauling2, Pauling3, Pauling4, Pauling5, RuleCannotBeAnalyzedError, \
    FrequencyEnvironmentPauling1, get_entropy_from_frequencies, get_most_frequent_environment
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from PlotClasses import PlotterPSE
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
import json
import scipy as sp
import numpy as np
from collections import OrderedDict
from collections import Counter
import matplotlib
import matplotlib.pyplot as plt
import os
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.analysis.chemenv.utils.scripts_utils import draw_cg
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.core.sites import PeriodicSite
from pymatgen.vis.structure_vtk import StructureVis
from yaml import load, dump
from pymatgen.core.periodic_table import Element

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

import csv


#TODO: discard non primitive structures, discard theoretical structures, update energy above hull for structures and discard structures not below or equal to 0.025 eV

class OverAllAnalysis:

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 analyse_structures=True, use_prematching=True, list_of_materials_to_investigate=None):
        """

        :param source:
        :param onlybinaries: only binary structures will be analysed
        :param plot_element_dependend_analysis:
        :param lowest_number_environments_for_plot: elements that have a lower number of environments in the evaluation will not be plotted
        :param lower_limit_plot:
        :param upper_limit_plot:
        :param analyse_structures: fulfiling and exceptional structures will be analysed
        :param use_prematching:
        :param list_of_materials_to_investigate:
        """

        self.source = source
        self.onlybinaries = onlybinaries
        self.plot_element_dependend_analysis = plot_element_dependend_analysis
        self.lowest_number_environments_for_plot = lowest_number_environments_for_plot
        self.analyse_structures = analyse_structures
        self.lower_limit_plot = lower_limit_plot
        self.upper_limit_plot = upper_limit_plot
        self.use_prematching = use_prematching
        self.list_of_materials_to_investigate = list_of_materials_to_investigate
        # could include a function that starts plotting from saved data?
        # how should one do this in the best way

    def _get_list_materials(self, source='MP', onlybinaries=False, start_material=None, stop_material=None):
        """

        :param source:
        :param onlybinaries: from the list of materials, only binaries are used
        :param start_material: cuts the list of all materials not the binary list
        :param stop_material: cuts the list of all materials not the binary list
        :return:
        """
        #TODO: test this part
        if source == 'MP':
            with open("../Auswertung/Should_not_be_changed/allmaterials.json", "r") as f:
                list_compound_dict = json.load(f)
                list_compound_dict["is_clear_compounds"] = list(set(list_compound_dict["is_clear_compounds"]))
                if not onlybinaries:
                    if start_material is None and stop_material is None:
                        list_compound = list_compound_dict["is_clear_compounds"]
                    elif start_material is None and stop_material is not None:
                        list_compound = list_compound_dict["is_clear_compounds"][0:stop_material]
                    elif start_material is not None and stop_material is None:
                        list_compound = list_compound_dict["is_clear_compounds"][start_material:]
                    else:
                        list_compound = list_compound_dict["is_clear_compounds"][start_material:stop_material]
                else:
                    list_compound=list_compound_dict["is_clear_compounds"]

        #TODO: generate new list of materials with these criteria!!!
        elif source == 'MP_very_symmetric':
            with open("../Auswertung/Should_not_be_changed/ce_fraction_0.95plus_csm_0.1plus_eabovehull0.025plus_discardS1.json",
                      "r") as f:
                list_compound_dict = json.load(f)
                if not onlybinaries:
                    if start_material is None and stop_material is None:
                        list_compound = list_compound_dict["is_clear_compounds"]
                    elif start_material is None and stop_material is not None:
                        list_compound = list_compound_dict["is_clear_compounds"][0:stop_material]
                    elif start_material is not None and stop_material is None:
                        list_compound = list_compound_dict["is_clear_compounds"][start_material:]
                    else:
                        list_compound = list_compound_dict["is_clear_compounds"][start_material:stop_material]
                else:
                    list_compound=list_compound_dict["is_clear_compounds"]
        #TODO: include test for start_material != None and stop_material==None
        elif source == 'my_own_list':
            list_compound = self._get_precomputed_results(self.list_of_materials_to_investigate)
            if not onlybinaries:
                if start_material is None and stop_material is None:
                    list_compound = list_compound
                elif start_material is None and stop_material is not None:
                    list_compound = list_compound[0:stop_material]
                elif start_material is not None and stop_material is None:
                    list_compound = list_compound[start_material:]
                else:
                    list_compound = list_compound[start_material:stop_material]
            else:
                list_compound = list_compound
        #TODO: test all possibilities
        # here: be careful that you are aware of start and stop_material !!
        if onlybinaries:
            list_compound_binaries = []
            for mat in list_compound:
                numberofcat = 0
                lse = self._get_lse_from_folder(mat=mat, source=source)
                for el in lse.structure.composition:
                    if str(el) != 'O':
                        numberofcat += 1
                if numberofcat == 1:
                    list_compound_binaries.append(mat)
            if start_material is None and stop_material is None:
                list_compound = list_compound_binaries
            elif start_material is None and stop_material is not None:
                list_compound = list_compound_binaries[0:stop_material]
            elif start_material is not None and stop_material is None:
                list_compound = list_compound_binaries[start_material:]
            else:
                list_compound = list_compound_binaries[start_material:stop_material]

            return list_compound

        else:
            return list_compound

    def _get_lse_from_folder(self, mat, source='MP'):
        if source == 'MP' or source == 'MP_very_symmetric' or source == 'my_own_list':
            with open("../../DB_chem_env/1st_rule/" + mat + ".json", 'r') as f:
                data = json.load(f)

        lse = LightStructureEnvironments.from_dict(data)
        return lse

    def _plot_PSE(self, Dict_to_Plot, lowest_number_of_environments_considered, xlim=[1, 18], ylim=[1, 10],
                  lowerlimit=0,
                  upperlimit=1, counter_cations_env=None, plot_directly_from_freq=False):


        plotterpse = PlotterPSE(valuestoplot=Dict_to_Plot, counter_cations_env=counter_cations_env,
                                plot_directly_from_freq=plot_directly_from_freq)
        # atoms to consider for plot instead?

        plt = plotterpse.get_plot(xlim=xlim, ylim=ylim,
                                  lowest_number_of_environments_considered=lowest_number_of_environments_considered,
                                  lowerlimit=lowerlimit, upperlimit=upperlimit)
        return plt

    def _save_results_to_file(self, dict_to_save, path):
        with open(path, 'w') as f:
            json.dump(dict_to_save, f)

    def _get_precomputed_results(self, path):
        with open(path, 'r') as f:
            dict_here = json.load(f)
        return dict_here

    def _add_dict_cat_dependency(self, start_dict, dict_to_add, number_of_elements_to_add=2):
        """

        :param start_dict: dict that will be extended
        :param dict_to_add: dict that will be used to extend
        :param number_of_elements_to_add:
             ==1: following dicts can be summed, integers for each key will be summed: {"Ga":1, "Sn":2}
             ==2: following dicts can be summed, individual numbers in the arrays of each key will be summed: {"Ga": [1,0], "Sn": [0,2]}
             ==3: does not exist at the moment
             ==4: following dicts can be extended with another dict: {"Ga":["O:6","O:6"], "Sn":["T:4","T:4"]}
        :return: None, start_dict will have the new values
        """

        for key, item in dict_to_add.items():
            if not key in start_dict:
                if number_of_elements_to_add == 1:
                    start_dict[key] = dict_to_add[key]
                elif number_of_elements_to_add == 2:
                    start_dict[key] = dict_to_add[key].copy()
                # elif number_of_elements_to_add == 3:
                #     start_dict[key] = dict_to_add[key].copy()
                elif number_of_elements_to_add == 4:
                    start_dict[key] = dict_to_add[key].copy()
            else:
                if number_of_elements_to_add == 1:
                    start_dict[key] += dict_to_add[key]
                elif number_of_elements_to_add == 2:
                    start_dict[key][0] += dict_to_add[key][0]
                    start_dict[key][1] += dict_to_add[key][1]
                # elif number_of_elements_to_add == 3:
                #     for env, number in dict_to_add[key].items():
                #         if env not in start_dict[key]:
                #             start_dict[key][env] = dict_to_add[key][env]
                #         else:
                #             start_dict[key][env] += dict_to_add[key][env]
                elif number_of_elements_to_add == 4:
                    start_dict[key].extend(dict_to_add[key])

    def _get_similar_structures(self, list_mat_id, source='MP', save_to_file=True,
                                path_to_save='Similar_Structures.json', fetch_results_only=False,
                                start_from_Matching=False,
                                path_to_precomputed_matching="Should_not_be_changed/Matching_All_Structures.json",
                                restart_from_matching=False):
        """

        :param list_mat_id:
        :param source:
        :param save_to_file:
        :param path_to_save:
        :param fetch_results_only:
        :param start_from_Matching: matches with the help of existing matching
        :param restart_from_matching: continues a matching from a file, e.g. to continue it later when interrupted
        :return:
        """

        if not start_from_Matching and not fetch_results_only:
            if not restart_from_matching:
                dictstructures = {}
                information_mat = {}
            else:
                outputdict = self._get_precomputed_results(path_to_save)
                dictstructures = outputdict['structure_matching']
                information_mat = outputdict['additional_info']

            for mat in list_mat_id:
                lse = self._get_lse_from_folder(mat, source=source)
                structure = lse.structure

                if mat not in information_mat:
                    information_mat[mat] = lse.structure.composition.reduced_formula
                if len(dictstructures.keys()) != 0:

                    Matcher = StructureMatcher(attempt_supercell=True, comparator=FrameworkComparator())

                    found = False
                    foundmat = ''
                    for mat2 in dictstructures:
                        lse2 = self._get_lse_from_folder(mat2, source=source)
                        structure2 = lse2.structure

                        result = Matcher.fit(structure, structure2)
                        if result:
                            found = True
                            foundmat = mat2
                            break;

                    if not found:
                        dictstructures[str(mat)] = []
                        dictstructures[str(mat)].append(mat)
                    else:
                        if mat not in dictstructures[str(foundmat)]:
                            dictstructures[str(foundmat)].append(mat)

                    dictstructures = OrderedDict(sorted(dictstructures.items(), key=lambda t: len(t[1]), reverse=True))

                else:
                    dictstructures[str(mat)] = []
                    dictstructures[str(mat)].append(mat)
            else:
                pass
        elif start_from_Matching:
            prematching = self._get_precomputed_results(path_to_precomputed_matching)

            information_mat = {}
            new_dict = {}
            for mat in list_mat_id:
                lse = self._get_lse_from_folder(mat, source=source)

                if mat not in information_mat:
                    information_mat[mat] = lse.structure.composition.reduced_formula

                in_there = False
                for key, values in prematching['structure_matching'].items():
                    if mat in values:
                        in_there = True
                        if not key in new_dict:
                            new_dict[key] = []
                        new_dict[key].append(mat)
                    # else:
                    #     raise ValueError
                if not in_there:
                    structure = lse.structure
                    Matcher = StructureMatcher(attempt_supercell=True, comparator=FrameworkComparator())
                    found = False
                    foundmat = ''
                    for mat2 in prematching['structure_matching'].keys():
                        lse2 = self._get_lse_from_folder(mat2, source=source)
                        structure2 = lse2.structure

                        result = Matcher.fit(structure, structure2)
                        if result:
                            found = True
                            foundmat = mat2
                            break;
                    if found:
                        #TODO: check if this can ever happen
                        if not foundmat in new_dict:
                            new_dict[foundmat] = []
                        if mat not in new_dict[foundmat]:
                            new_dict[foundmat].append(mat)
                    else:
                        if not mat in new_dict:
                            new_dict[mat] = []
                        if mat not in new_dict[mat]:
                            new_dict[mat].append(mat)
            dictstructures = {}
            for values in new_dict.values():
                dictstructures[values[0]] = values

            dictstructures = OrderedDict(sorted(dictstructures.items(), key=lambda t: len(t[1]), reverse=True))

        if not fetch_results_only:
            new_list_mat_id = []
            for key, item in dictstructures.items():
                for item2 in item:
                    if item2 not in new_list_mat_id:
                        new_list_mat_id.append(item2)

            outputdict = {}
            outputdict['list_mat_id'] = new_list_mat_id
            outputdict['structure_matching'] = dictstructures
            outputdict['additional_info'] = information_mat

        if save_to_file and (not fetch_results_only) and (path_to_save is not None):
            self._save_results_to_file(outputdict, path_to_save)
        #TODO: test valueerror -> use another list before loading (start and stop_material nutzen)
        if fetch_results_only:
            outputdict = self._get_precomputed_results(path_to_save)
            if not set(outputdict['list_mat_id']) == set(list_mat_id):
                raise ValueError
        return outputdict

    # def _print_to_screen_similar_structures(self, dict_to_print, add_properties_from_lse=False, source='MP'):
    #     for key, items in OrderedDict(
    #             sorted(dict_to_print['structure_matching'].items(), key=lambda t: len(t[1]), reverse=True)).items():
    #
    #         if len(items) == 1:
    #             print(key + ' (' + str(dict_to_print['additional_info'][key]) + ', ' + str(
    #                 len(items)) + ' representative): ', end='')
    #         else:
    #             print(key + ' (' + str(dict_to_print['additional_info'][key]) + ', ' + str(
    #                 len(items)) + ' representatives): ', end='')
    #         for item in items:
    #             if add_properties_from_lse:
    #                 lse = self._get_lse_from_folder(mat=item, source=source)
    #                 valence = [i for i in lse.valences if i > 0]
    #                 # get CN:
    #                 CN = []
    #                 for isite, site_envs in enumerate(lse.coordination_environments):
    #
    #                     if site_envs != None:
    #                         if len(site_envs) > 0:
    #                             CN.append(site_envs[0]['ce_symbol'].split(':')[1])
    #                 print(item + ' (' + str(dict_to_print['additional_info'][item]) + ', CNs:' + str(
    #                     CN) + ', valences: ' + str(valence) + '), ', end='')
    #             print(item + ' (' + str(dict_to_print['additional_info'][item]) + '), ', end='')
    #
    #         print()

    def _print_to_file_similar_structures(self, dict_to_print, source='MP', filename='Test.yaml', fmt='yml'):
        if fmt == 'yml':
            new_dict_to_print = OrderedDict()
            for key, items in OrderedDict(
                    sorted(dict_to_print['structure_matching'].items(), key=lambda t: len(t[1]), reverse=True)).items():

                new_dict_to_print[key] = {}

                for item in items:

                    lse = self._get_lse_from_folder(mat=item, source=source)
                    valence = [i for i in lse.valences if i > 0]
                    # get CN:
                    cat = []
                    CN = []
                    for isite, site_envs in enumerate(lse.coordination_environments):

                        if site_envs != None:
                            if len(site_envs) > 0:
                                CN.append(site_envs[0]['ce_symbol'].split(':')[1])
                                cat.append(lse.structure[isite].species_string)
                    new_dict_to_print[key][item] = {"formula": dict_to_print['additional_info'][item], "CN": CN,
                                                    "valences": valence, "cations": cat}

            with open(filename, 'w') as f:
                dump(OrderedDict(new_dict_to_print), f)

        elif fmt == 'csv':
            new_dict_to_print = []
            for key, items in OrderedDict(
                    sorted(dict_to_print['structure_matching'].items(), key=lambda t: len(t[1]), reverse=True)).items():

                for item in items:

                    lse = self._get_lse_from_folder(mat=item, source=source)
                    valence = [i for i in lse.valences if i > 0]
                    # get CN:
                    cat = []
                    CN = []
                    for isite, site_envs in enumerate(lse.coordination_environments):

                        if site_envs != None:
                            if len(site_envs) > 0:
                                CN.append(site_envs[0]['ce_symbol'].split(':')[1])
                                cat.append(lse.structure[isite].species_string)
                    new_dict_to_print.append(
                        {"mpid": str(item), "formula": dict_to_print['additional_info'][item], "CN": CN,
                         "valences": valence, "cations": cat, "structure_type": key})

            with open(filename, 'w') as csvfile:
                filewriter = csv.writer(csvfile, delimiter=';',
                                        quotechar='|', quoting=csv.QUOTE_MINIMAL)
                filewriter.writerow(['mp-id', 'formula', 'structure_type', 'cations', 'valences', 'CNs'])
                for line in new_dict_to_print:
                    filewriter.writerow(
                        [str(line["mpid"]), str(line["formula"]), str(line["structure_type"]), str(line["cations"]),
                         str(line["valences"]), str(line["CN"])])

    def _pieplot_connections(self, corner, edge, face, title="All"):

        labels = 'via corners', 'via edges', 'via faces'
        sizes = [corner, edge, face]
        colors = ['#5da5daff', '#faa43aff', '#f15854ff']

        plt.pie(sizes, labels=labels, colors=colors,
                autopct='%1.1f%%', shadow=False, startangle=140)
        plt.title(title)
        plt.axis('equal')
        return plt

    def visualize_structure_by_id(self, mat, supercell=[1, 1, 0]):
        # TODO: implement supercell
        # get information on the supercell to neigbors

        lse = self._get_lse_from_folder(mat, source=self.source)

        vis = StructureVis(show_polyhedron=False, show_unit_cell=True)
        vis.show_help = False
        vis.set_structure(lse.structure)

        # look how to implement the supercell correctly
        # strategy = lse.strategy
        allcg = AllCoordinationGeometries()

        for isite, site in enumerate(lse.structure):
            if lse.coordination_environments[isite] is None:
                continue
            try:
                ces = lse.coordination_environments[isite]
            except NeighborsNotComputedChemenvError:
                continue
            if len(ces) == 0:
                continue
            ce = ces[0]
            # print(ce)
            newdeltas1 = [[float(i), float(0), float(0)] for i in range(1, supercell[0] + 1)]
            newdeltas2 = [[float(0), float(i), float(0)] for i in range(1, supercell[1] + 1)]
            newdeltas3 = [[float(0), float(0), float(i)] for i in range(1, supercell[2] + 1)]
            newdeltas = newdeltas1 + newdeltas2 + newdeltas3 + [0.0, 0.0, 0.0]
            print(newdeltas)

            mydeltas = newdeltas
            for mydelta in mydeltas:
                psite = PeriodicSite(site.species, site.frac_coords + mydelta, site.lattice,
                                     properties=site.properties)
                vis.add_site(psite)
                neighbors_start = lse.neighbors_sets[isite][0].neighb_sites
                neighbors = []
                for neighbor in neighbors_start:
                    neighbors.append(PeriodicSite(neighbor.species, neighbor.frac_coords + mydelta, neighbor.lattice,
                                                  properties=neighbor.properties))


                draw_cg(vis, psite, neighbors, cg=allcg.get_geometry_from_mp_symbol(ce['ce_symbol']),
                        perm=ce['permutation'])
        vis.show()


# Not sure this should be tested - let's see
# class MatchAllStructures(OverAllAnalysis):
#     def __init__(self, source='MP', startmaterial=None, stopmaterial=None, restart_from_matching=False):
#         self.source = source
#         materials_list = self._get_list_materials(source=self.source, start_material=startmaterial,
#                                                   stop_material=stopmaterial)
#         self._get_similar_structures(materials_list, source=source, save_to_file=True,
#                                      path_to_save="Results/Matching_All_Structures.json",
#                                      restart_from_matching=restart_from_matching)
#
#
# # most frequent environment: this is also something that needs to be computed
# might need to be included in PaulingRules


class Pauling1Frequency(OverAllAnalysis):
    """Analyses the 1-percentropy of certain cation"""

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 list_of_materials_to_investigate=None,start_material=None,stop_material=None):
        """
            :param source:
            :param onlybinaries: only binary structures will be analysed
            :param analyse_structures: fulfiling and exceptional structures will be analysed
            :param lowest_number_environments_for_plot: elements that have a lower number of environments in the evaluation will not be plotted
            :param save_result_data: all results will be saved so that results can easily be recreated
            :return:
        """

        self.source = source
        self.onlybinaries = onlybinaries
        self.plot_element_dependend_analysis = plot_element_dependend_analysis
        self.lowest_number_environments_for_plot = lowest_number_environments_for_plot
        self.lower_limit_plot = lower_limit_plot
        self.upper_limit_plot = upper_limit_plot
        self.start_material=start_material
        self.stop_material=stop_material
        self.list_of_materials_to_investigate = list_of_materials_to_investigate

    def run(self, start_from_results=False, save_result_data=True, path_to_save='Results/Results_First_Limits.json'):

        if not start_from_results:
            self._new_setup()
            self.Plot_PSE_most_frequent = get_most_frequent_environment(self.All_Details)


        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.Plot_PSE_most_frequent = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])

        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_most_frequent,
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 plot_directly_from_freq=True, lowerlimit=self.lower_limit_plot,
                                 upperlimit=self.upper_limit_plot, counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['PSE_Dict'] = self.Plot_PSE_most_frequent
            outputdict['Counter_cation'] = dict(self.present_env)
            self._save_results_to_file(outputdict, path_to_save)

    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,start_material=self.start_material,stop_material=self.stop_material)
        self.All_Details = {}
        self.present_env = {}
        for mat in list_mat:
            # print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)
            pauling0 = Pauling0(lse)
            pauling1_limit = FrequencyEnvironmentPauling1(lse=lse)
            Details = pauling1_limit.get_details()

            self._add_dict_cat_dependency(start_dict=self.All_Details, dict_to_add=Details, number_of_elements_to_add=4)
            self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                          number_of_elements_to_add=1)


class Pauling1Entropy(Pauling1Frequency):
    """Analyses the 1-percentropy of certain cation"""

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 list_of_materials_to_investigate=None,start_material=None,stop_material=None):
        """
            :param source:
            :param onlybinaries: only binary structures will be analysed
            :param analyse_structures: fulfiling and exceptional structures will be analysed
            :param lowest_number_environments_for_plot: elements that have a lower number of environments in the evaluation will not be plotted
            :param save_result_data: all results will be saved so that results can easily be recreated
            :return:
        """

        super(Pauling1Entropy, self).__init__(source=source, onlybinaries=onlybinaries,
                                              plot_element_dependend_analysis=plot_element_dependend_analysis,
                                              lowest_number_environments_for_plot=lowest_number_environments_for_plot,
                                              lower_limit_plot=lower_limit_plot, upper_limit_plot=upper_limit_plot,
                                              list_of_materials_to_investigate=list_of_materials_to_investigate,start_material=start_material,stop_material=stop_material)

    #     self.source = source
    #     self.onlybinaries = onlybinaries
    #     self.plot_element_dependend_analysis = plot_element_dependend_analysis
    #     self.lowest_number_environments_for_plot = lowest_number_environments_for_plot
    #     self.lower_limit_plot = lower_limit_plot
    #     self.upper_limit_plot = upper_limit_plot
    #     if list_of_materials_to_investigate is not None:
    #         self.list_of_materials_to_investigate=list_of_materials_to_investigate

    def run(self, start_from_results=False, save_result_data=True, path_to_save='Results/Results_First_Limits.json',start_material=None,stop_material=None):


        if not start_from_results:
            super(Pauling1Entropy, self)._new_setup()
            self.Plot_PSE_entropy = get_entropy_from_frequencies(self.All_Details)

        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.Plot_PSE_entropy = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])

        if self.plot_element_dependend_analysis:

            plt = self._plot_PSE(self.Plot_PSE_entropy,
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 plot_directly_from_freq=True, lowerlimit=self.lower_limit_plot,
                                 upperlimit=self.upper_limit_plot, counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['PSE_Dict'] = self.Plot_PSE_entropy
            outputdict['Counter_cation'] = dict(self.present_env)
            self._save_results_to_file(outputdict, path_to_save)

    # def _new_setup(self):
    #     list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,start_material=self.start_material,stop_material=self.stop_material)
    #     self.All_Details = {}
    #     self.present_env = {}
    #     for mat in list_mat:
    #         # print(mat)
    #         lse = self._get_lse_from_folder(mat, source=self.source)
    #         pauling0 = Pauling0(lse)
    #         pauling1_limit = FrequencyEnvironmentPauling1(lse=lse)
    #         Details = pauling1_limit.get_details()
    #
    #         self._add_dict_cat_dependency(start_dict=self.All_Details, dict_to_add=Details, number_of_elements_to_add=4)
    #         self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
    #                                       number_of_elements_to_add=1)


class Pauling1OverAllAnalysis(OverAllAnalysis):

    def run(self, start_from_results=False, save_result_data=True, save_structure_analysis=True,
            restart_from_saved_structure_analyisis=False,
            path_to_save='Results/Results_First_Rule.json', start_material=None,stop_material=None):

        if not start_from_results:
            self.start_material = start_material
            self.stop_material = stop_material
            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']

        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_DICT,
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['PSE_Dict'] = self.Plot_PSE_DICT
            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split(
                                                                                 '.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)
            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split(
                                                                                 '.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling

            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _new_setup(self):
        # could include: newsetup or from file here!
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.Plot_PSE_DICT = {}
        # valence dependency can be introduced later

        for mat in list_mat:
            # print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)

            pauling1 = Pauling1(lse=lse, filenameradii='../univalent_cat_radii.json', onlylowerlimit=False)
            try:
                if pauling1.is_fulfilled():
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)
                # print('Cannot be analyzed')

            # get details and add them
            Details = pauling1.get_details()['cat_dependency']
            self._add_dict_cat_dependency(self.Plot_PSE_DICT, Details)


class Pauling2OverAllAnalysis(OverAllAnalysis):

    def run(self, show_plot=True, start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analyisis=False,
            save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json', start_material=None,stop_material=None):

        if not start_from_results:
            self.start_material = start_material
            self.stop_material = stop_material
            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])
            self.arraydev_share = np.array(inputdict['arraydev_share'])
            self.relativefrequency = np.array(inputdict['relativefrequency'])
            self.bs_sum_mean = inputdict['Mean_bvs']
            self.tot_stddev = np.array(inputdict['rel_std_dev_bvs'])

        if show_plot:
            plot = self._secondrule_plot(arraydev=self.arraydev_share * 100.0,
                                         relativefreqarray=self.relativefrequency * 100.0,
                                         tot_stddev=self.tot_stddev * 100.0)
            plot.show()

        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_DICT, xlim=[1, 18], ylim=[1, 10],
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 lowerlimit=self.lower_limit_plot, upperlimit=self.upper_limit_plot,
                                 counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['PSE_Dict'] = self.Plot_PSE_DICT
            outputdict['Counter_cation'] = dict(self.present_env)
            outputdict['arraydev_share'] = list(self.arraydev_share)
            outputdict['relativefrequency'] = list(self.relativefrequency)
            outputdict['Mean_bvs'] = self.bs_sum_mean
            outputdict['rel_std_dev_bvs'] = list(self.tot_stddev)
            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             start_from_Matching=self.use_prematching,
                                                                             fetch_results_only=restart_from_saved_structure_analyisis)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             start_from_Matching=self.use_prematching,
                                                                             fetch_results_only=restart_from_saved_structure_analyisis)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling

            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _stddev(self, lst, mean):
        """percental standard deviation of a sample in a list format"""
        # TODO: get rid of the list format
        sum = 0.0
        mn = mean
        for i in range(len(lst)):
            sum += pow((lst[i] - mn), 2)
        return np.sqrt([sum / (len(lst) - 1)]) / mn

    # def _stddev_normal(self, lst, mean):
    #     """?"""
    #     """standard deviation of a sample"""
    #     sum = 0.0
    #     mn = mean
    #     for i in range(len(lst)):
    #         sum += pow((lst[i] - mn), 2)
    #     return np.sqrt([sum / (len(lst) - 1)])

    def _get_deviation_from_ideal_value(self, inputarray, ideal):
        return [abs(float(x) - float(ideal)) for x in inputarray]

    def _get_frequency_of_values(self, arraydeviations, stepsize):
        frequency = []
        dev_array = np.arange(0, max(arraydeviations) + stepsize, stepsize)

        for step in dev_array:
            # print(step)

            frequency.append(len([x for x in arraydeviations if x <= step]))
            # print(frequency)

        return dev_array, frequency

    def _list_to_np_array_and_divide_by_value(self, inputarray, valuetodivide):
        return (np.array(inputarray) / float(valuetodivide))

    def _secondrule_plot(self, arraydev, relativefreqarray, tot_stddev, maxpercentage=70, save_plot=False,
                         filename='second_rule.svg'):
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        font = {'size': 22}

        matplotlib.rc('font', **font)
        plt.plot((arraydev), relativefreqarray)
        plt.plot([tot_stddev, tot_stddev], [0, 100.0], linewidth=3.0)
        plt.ylim([0, 100.0])
        plt.xlim([0, maxpercentage])
        plt.xticks(range(0, maxpercentage + 10, 10))
        plt.yticks(range(0, 100 + 20, 20))
        plt.xlabel("Absolute Deviation (%) from the ideal valence -2")
        plt.ylabel("Oxygen Atoms (%)")
        if save_plot:
            plt.savefig(filename)
        return plt

    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,start_material=self.start_material,stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        array_bvs = []
        self.Plot_PSE_DICT = {}
        self.present_env = {}
        # valence dependency can be introduced later
        for mat in list_mat:
            # print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)

            pauling2 = Pauling2(lse=lse)

            # one could also put this in a method!
            if pauling2.is_fulfilled():
                self.structures_fulfillingrule.append(mat)
            else:
                self.structures_exceptions.append(mat)

            # get details and add them
            Details = pauling2.get_details()
            bvs = Details['bvs_for_each_anion']

            self._add_dict_cat_dependency(self.Plot_PSE_DICT, Details['elementwise_fulfillment'])

            # add cations to each other, to count total number of cations
            self._add_dict_cat_dependency(self.present_env, Details['cations_in_structure'],
                                          number_of_elements_to_add=1)

            array_bvs.extend(bvs)

            # print(self.array_bvs)

        self.array_bvs = array_bvs
        self.bs_sum_mean = np.mean(array_bvs)

        self.tot_stddev = self._stddev(array_bvs, np.mean(array_bvs))
        ideal_bs = 2.0
        bs_dev = self._get_deviation_from_ideal_value(inputarray=array_bvs, ideal=ideal_bs)
        self.bs_dev = bs_dev
        dev_array, self.frequency = self._get_frequency_of_values(arraydeviations=bs_dev, stepsize=0.01)

        self.arraydev_share = self._list_to_np_array_and_divide_by_value(np.array(dev_array), ideal_bs)
        self.relativefrequency = self._list_to_np_array_and_divide_by_value(np.array(self.frequency),
                                                                            float(len(array_bvs)))



class Pauling3OverAllAnalysis(OverAllAnalysis):

    def run(self, show_plot=True, start_from_connections=False, save_connections=True,
            connections_folder='AnalysisConnections', start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analyisis=False,  save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json',start_material=None,stop_material=None):

        self.start_material = start_material
        self.stop_material = stop_material
        self.connections_folder = connections_folder
        self.start_from_connections = start_from_connections
        self.save_connections = save_connections


        if not start_from_results:
            self._new_setup()
        else:
            pass
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])
            self.connections = inputdict['connections']

        if show_plot:
            self._plot_influence_atomic_radii(self.Plot_PSE_DICT)

            exit()
            plot = self._pieplot_connections(self.connections['corner'], self.connections['edge'],
                                             self.connections['face'], 'Connected Pairs of Polyhedra')
            plot.show()


            #TODO: plot that shows number of connections via corners and edges vs. atomic radii, and number of connections via faces vs. atomic radii: get radii first
        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_DICT, xlim=[1, 18], ylim=[1, 10],
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 lowerlimit=self.lower_limit_plot, upperlimit=self.upper_limit_plot,
                                 counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['PSE_Dict'] = self.Plot_PSE_DICT
            outputdict['Counter_cation'] = dict(self.present_env)
            outputdict['connections'] = self.connections

            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling


            if save_structure_analysis:
                #self._print_to_screen_similar_structures(dict_similarstructures_fulfilling)
                self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                       0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                       0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _plot_influence_atomic_radii(self,Plot_PSE_Dict):
        #oder histogramm
        corneredge_histo=[]
        face_histo=[]
        for key,item in Plot_PSE_Dict.items():
            for number in range(0,item[0]):
                corneredge_histo.append(Element(key).atomic_radius)
            for number in range(0,item[1]):
                face_histo.append(Element(key).atomic_radius)

        print(np.mean(corneredge_histo))
        print(np.std(corneredge_histo,ddof=1))
        print(np.mean(face_histo))
        print(np.std(face_histo,ddof=1))

        #range noch korrekt anpassen
        plt.subplot(2,1,1)
        plt.hist(corneredge_histo,bins=len(range(int(min(corneredge_histo)*100),int(max(corneredge_histo)*100)+5)),align='mid',alpha=1)
        plt.axvline(np.mean(corneredge_histo),color='r')

        #plt.plot([0, 1000], [np.mean(corneredge_histo), np.mean(corneredge_histo)], 'b-')
        # n,bins,patches=plt.hist(x=corneredge_histo,bins=len(range(int(min(corneredge_histo)*100),int(max(corneredge_histo)*100)+5)))
        # y =sp.stats.norm.pdf(bins,np.mean(corneredge_histo),np.std(corneredge_histo,ddof=1))
        # plt.plot(bins,y,'r--',linewidth=1)
        #
        plt.subplot(2,1,2)
        plt.hist(face_histo, bins=len(range(int(min(face_histo) * 100), int(max(face_histo) * 100) + 5)), align='mid',
                 alpha=1)
        plt.axvline(np.mean(face_histo), color='r')
        #plt.plot([0,1000],[np.mean(face_histo),np.mean(face_histo)],'b-')
        plt.show()
    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.Plot_PSE_DICT = {}
        self.connections = {"corner": 0, "edge": 0, "face": 0}
        self.present_env = {}
        # valence dependency can be introduced later
        #TODO: should search for an exception to third rule
        for mat in list_mat:
            #print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)
            pauling0 = Pauling0(lse)
            pauling3 = Pauling3()
            if not self.start_from_connections or not os.path.isfile(
                    os.path.join(self.connections_folder, mat + '.json')):
                pauling3.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder, distance=8.0)
            else:
                pauling3.from_file(filename=mat + '.json', foldername=self.connections_folder)
            try:
                if pauling3.is_fulfilled():
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)

            # get details and add them

            # weiteres dict anlegen, um verknuepfungsplot zu machen

            Details = pauling3.get_details()
            # print(Details)
            New_Details = self._reformat_details(Details['species'])
            self._add_dict_cat_dependency(self.Plot_PSE_DICT, New_Details)

            self._sum_connections(self.connections, Details)

            self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                          number_of_elements_to_add=1)

            # function that brings Details in correct form:
            # for all keys, sum over all valences, add 'corner'+'edge' and all 'face together,
            # give an output dict that can be used for PlotPSE

            # bvs = Details['bvs_for_each_anion']

        # print(self.Plot_PSE_DICT)
        # print(self.connections)
        # print(self.present_env)
        # add cations to each other, to count total number of cations
        # self._add_dict_cat_dependency(self.present_env, Details['cations_in_structure'],
        #                              number_of_elements_to_add=1)

        # array_bvs.extend(bvs)

        # print(self.array_bvs)

    def _reformat_details(self, Details):
        New_Details = {}
        for key, item in Details.items():
            if not key in New_Details:
                New_Details[key] = [0, 0]
            for item2 in item.values():
                New_Details[key][0] += item2['corner']
                New_Details[key][0] += item2['edge']
                New_Details[key][1] += item2['face']
        return New_Details

    def _sum_connections(self, start_dict, Details):
        start_dict["corner"] += Details["corner"]
        start_dict["edge"] += Details["edge"]
        start_dict["face"] += Details["face"]

        pass


class Pauling4OverAllAnalysis(OverAllAnalysis):
    #TODO: think about another elementwise analysis instead -> could also just count for fulfilling or unfulfilling for cations that have lowest CN and highest valence

    def run(self, show_plot=True, start_from_connections=False, save_connections=True,
            connections_folder='AnalysisConnections', start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analyisis=False,save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json', start_material=None, stop_material=None):

        self.connections_folder = connections_folder
        self.start_from_connections = start_from_connections
        self.save_connections = save_connections

        # for generating connections files only:
        self.start_material = start_material
        self.stop_material = stop_material

        if not start_from_results:
            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])
            self.Dict_val = inputdict['val_dep_for_plot']
            self.Dict_CN = inputdict['CN_dep_for_plot']

        if show_plot:

            plot = self._fourthrule_plot(self.Dict_val, option='val', vmin=0.70, vmax=1.0)
            plot.show()
            plot = self._fourthrule_plot(self.Dict_CN, option='CN', vmin=0.4, vmax=1.0)
            plot.show()


        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_DICT, xlim=[1, 18], ylim=[1, 10],
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 lowerlimit=self.lower_limit_plot, upperlimit=self.upper_limit_plot,
                                 counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['PSE_Dict'] = self.Plot_PSE_DICT
            outputdict['Counter_cation'] = dict(self.present_env)
            outputdict['val_dep_for_plot'] = self.Dict_val
            outputdict['CN_dep_for_plot'] = self.Dict_CN
            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling


            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                       0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                       0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.Dict_val = {}
        self.Dict_CN = {}
        self.Plot_PSE_DICT = {}
        self.present_env = {}
        # valence dependency can be introduced later

        for mat in list_mat:
            lse = self._get_lse_from_folder(mat, source=self.source)

            pauling4 = Pauling4()
            if not self.start_from_connections or not os.path.isfile(
                    os.path.join(self.connections_folder, mat + '.json')):
                pauling4.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder, distance=8.0)
            else:
                pauling4.from_file(filename=mat + '.json', foldername=self.connections_folder)
            try:
                if pauling4.is_fulfilled():
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
                pauling0 = Pauling0(lse)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)

            try:
                Details = pauling4.get_details()


                # think about symmetry - there might be something wrong here
                New_Details_val = self._reformat_details_val(Details)
                New_Details_CN = self._reformat_details_CN(Details)

                # # tests if matrices are symmetric!
                # if not self._check_details_symmetric(New_Details_val):
                #     raise ValueError
                # if not self._check_details_symmetric(New_Details_CN):
                #     raise ValueError

                self._add_dict_CN_val(self.Dict_val, New_Details_val)
                self._add_dict_CN_val(self.Dict_CN, New_Details_CN)

                Elementwise_analysis = self._reformat_details_elementwise(Details)

                self._add_dict_cat_dependency(self.Plot_PSE_DICT, Elementwise_analysis)
                #print(Details['elementwise'])
                #print(Details['maxval'])
                #print(Details['minCN'])
                #print(Elementwise_analysis)

                self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                              number_of_elements_to_add=1)

            except RuleCannotBeAnalyzedError:
                pass

    def _reformat_details_elementwise(self, Details):
        maxval1 = 'val1:' + str(Details['maxval'])
        maxval2 = 'val2:' + str(Details['maxval'])
        minCN1 = 'CN1:' + str(Details['minCN'])
        minCN2 = 'CN2:' + str(Details['minCN'])

        outputdict = {}
        for el, item in Details['elementwise'].items():
            if maxval1 in item:
                if maxval2 in item[maxval1]:
                    if minCN1 in item[maxval1][maxval2]:
                        if minCN2 in item[maxval1][maxval2][minCN1]:
                            connections = item[maxval1][maxval2][minCN1][minCN2]['corner'] + \
                                          item[maxval1][maxval2][minCN1][minCN2]['edge'] + \
                                          item[maxval1][maxval2][minCN1][minCN2]['face']
                            if connections > 0:
                                if el not in outputdict:
                                    outputdict[el] = [0, 0]
                                outputdict[el][0] += item[maxval1][maxval2][minCN1][minCN2]['no']
                                outputdict[el][1] += connections
                            else:
                                if el not in outputdict:
                                    outputdict[el] = [0, 0]
                                outputdict[el][0] += item[maxval1][maxval2][minCN1][minCN2]['no']
                                outputdict[el][1] += connections

        return outputdict

    def _add_dict_CN_val(self, start_dict, to_add_dict):
        for key, item in to_add_dict.items():
            if key not in start_dict:
                start_dict[key] = {}
            for key2, item2 in item.items():
                if key2 not in start_dict[key]:
                    start_dict[key][key2] = [0, 0]
                start_dict[key][key2][0] += to_add_dict[key][key2][0]
                start_dict[key][key2][1] += to_add_dict[key][key2][1]

    def _reformat_details_val(self, Details):
        New_Details = {}

        for key1, item1 in Details.items():

            if "val1" in key1:
                val1 = key1.split(":")[1]
                if val1 not in New_Details:
                    New_Details[val1] = {}
                for key2, item2 in item1.items():
                    val2 = key2.split(":")[1]
                    if val2 not in New_Details[val1]:
                        New_Details[val1][val2] = [0, 0]
                    if val2 not in New_Details:
                        New_Details[val2] = {}
                    if val1 not in New_Details[val2]:
                        New_Details[val2][val1] = [0, 0]
                    for key3, item3 in item2.items():
                        for key4, item4 in item3.items():
                            New_Details[val1][val2][0] += item4["no"]
                            New_Details[val1][val2][1] += item4["corner"]
                            New_Details[val1][val2][1] += item4["edge"]
                            New_Details[val1][val2][1] += item4["face"]
                            if val1!=val2:
                                New_Details[val2][val1][0] += item4["no"]
                                New_Details[val2][val1][1] += item4["corner"]
                                New_Details[val2][val1][1] += item4["edge"]
                                New_Details[val2][val1][1] += item4["face"]

        return New_Details

    # def _check_details_symmetric(self, inputdict):
    #     # teste, ob ausgabe symmetrisch
    #     for key1, item in inputdict.items():
    #         for key2, item2 in item.items():
    #             if key1 <= key2:
    #                 if not item2 == inputdict[key2][key1]:
    #                     return False
    #     return True

    def _reformat_details_CN(self, Details):
        New_Details = {}

        for key1, item1 in Details.items():

            if "val1" in key1:
                # val1 = key1.split(":")[1]
                for key2, item2 in item1.items():
                    # val2 = key2.split(":")[1]
                    for key3, item3 in item2.items():
                        CN1 = key3.split(":")[1]
                        if CN1 not in New_Details:
                            New_Details[CN1] = {}

                        for key4, item4 in item3.items():
                            CN2 = key4.split(":")[1]
                            if CN2 not in New_Details[CN1]:
                                New_Details[CN1][CN2] = [0, 0]
                            if CN2 not in New_Details:
                                New_Details[CN2] = {}
                            if CN1 not in New_Details[CN2]:
                                New_Details[CN2][CN1] = [0, 0]
                            New_Details[CN1][CN2][0] += item4["no"]
                            New_Details[CN1][CN2][1] += item4["corner"]
                            New_Details[CN1][CN2][1] += item4["edge"]
                            New_Details[CN1][CN2][1] += item4["face"]
                            if CN1!=CN2:
                                New_Details[CN2][CN1][0] += item4["no"]
                                New_Details[CN2][CN1][1] += item4["corner"]
                                New_Details[CN2][CN1][1] += item4["edge"]
                                New_Details[CN2][CN1][1] += item4["face"]

        return New_Details

    def _plot_get_max(self,inputdict):
        return  max([int(x) for x in inputdict.keys()])

    def _plot_init_np(self,maxvalue):
        return np.full((maxvalue,maxvalue),np.nan)

    def _plot_fill_np(self,PlotEndDict,inputdict,lowervalue):
        for key, item in inputdict.items():
            for key2, item2 in item.items():
                if item2[0] + item2[1] > lowervalue:
                    PlotEndDict[int(key)][int(key2)] = float(item2[0]) / float(item2[0] + item2[1])


    def _fourthrule_plot(self, inputdict, option='CN', vmin=0.0, vmax=1.0, lowervalue=50, leaveoutCN13=True):
        """
        #leaves out CN 13 for option "CN" -> too few entries
        :param inputdict:
        :param option: 'CN' or 'val' -> changes labels of axis dynamically
        :return:
        """
        # get valmax
        valCNmax = self._plot_get_max(inputdict)

        if option == 'val':
            PlotEndDict = self._plot_init_np(valCNmax+1)
        elif option == 'CN':
            if not leaveoutCN13:
                PlotEndDict = self._plot_init_np(valCNmax+2)
            else:
                PlotEndDict = self._plot_init_np(valCNmax+2)

        self._plot_fill_np(PlotEndDict,inputdict,lowervalue)

        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        matplotlib.rcParams['figure.dpi'] = 800
        matplotlib.rc("savefig", dpi=800)
        font = {'family': 'normal',
                'size': 3}

        matplotlib.rc('font', **font)

        fig, (ax1) = plt.subplots(figsize=(1, 1), ncols=1)

        cmap = plt.cm.get_cmap('cool', 48)
        cmap.set_bad(color='white')
        cmap.set_under(color='black')

        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)

        pos = ax1.imshow(PlotEndDict, cmap=cmap, norm=norm)

        if option == 'CN':
            ax1.set_xlabel('CN1')

            if not leaveoutCN13:
                ax1.set_xticks(range(2, valCNmax + 2))
                ax1.set_xticklabels(range(2, valCNmax + 2))
                ax1.set_xlim(1.5, valCNmax + 1 - 0.5)
                ax1.set_yticks(range(2, valCNmax + 2))
                ax1.set_ylim(1.5, valCNmax + 1 - 0.5)
                ax1.set_yticklabels(range(2, valCNmax + 2))
            else:
                ax1.set_xticks(range(2, valCNmax + 1))
                ax1.set_xticklabels(range(2, valCNmax + 1))
                ax1.set_xlim(1.5, valCNmax + - 0.5)
                ax1.set_yticks(range(2, valCNmax + 1))
                ax1.set_ylim(1.5, valCNmax + - 0.5)
                ax1.set_yticklabels(range(2, valCNmax + 1))
            ax1.set_ylabel('CN')
        elif option == 'val':
            ax1.set_xlabel('val1')
            ax1.set_xticks(range(1, valCNmax + 1))
            ax1.set_xticklabels(range(1, valCNmax + 1))
            ax1.set_xlim(0.5, valCNmax - 0.5)

            ax1.set_yticks(range(1, valCNmax + 1))
            ax1.set_ylim(0.5, valCNmax - 0.5)
            ax1.set_yticklabels(range(1, valCNmax + 1))
            ax1.set_ylabel('val2')

        cb = fig.colorbar(pos, ax=ax1, cmap=cmap,
                          norm=norm,
                          orientation='vertical')
        ax1.set_title('Not-connected polyhedra')
        return plt


class Pauling5OverAllAnalysis(OverAllAnalysis):

    def run(self, remove_elements_low_entropy=True, show_plot=True, start_from_connections=False, save_connections=True,
            connections_folder='AnalysisConnections_5thRule', start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analyisis=False, save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json',
            threshold_remove_elements=0.95,start_material=None,stop_material=None):
        # threshold_remove_elements: 1-entropy that removes elements with only very few environments

        self.remove_elements_low_entropy = remove_elements_low_entropy
        self.connections_folder = connections_folder
        self.start_from_connections = start_from_connections
        self.save_connections = save_connections
        self.start_material = start_material
        self.stop_material = stop_material


        if not start_from_results:
            if self.remove_elements_low_entropy:
                newclass = Pauling1Entropy(source=self.source, onlybinaries=self.onlybinaries,
                                           plot_element_dependend_analysis=False, list_of_materials_to_investigate=self.list_of_materials_to_investigate)
                newclass.run(start_from_results=False, save_result_data=False)
                self.list_to_remove = [el for el, item in newclass.Plot_PSE_entropy.items() if
                                       item > threshold_remove_elements]
            else:
                self.list_to_remove=[]
            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])
            self.list_to_remove=inputdict['elements_low_entropy']

        if show_plot:
            plot = self._fifth_rule_plot(okay=len(self.structures_fulfillingrule),
                                         notokay=len(self.structures_exceptions))
            plot.show()

        if self.remove_elements_low_entropy and not start_from_results:

            new_env = {}
            for key, item in self.present_env.items():
                if key not in self.list_to_remove:
                    new_env[key] = item

        if self.plot_element_dependend_analysis:

            plt = self._plot_PSE(self.Plot_PSE_DICT, xlim=[1, 18], ylim=[1, 10],
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 lowerlimit=self.lower_limit_plot, upperlimit=self.upper_limit_plot,
                                 counter_cations_env=new_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['PSE_Dict'] = self.Plot_PSE_DICT
            outputdict['Counter_cation'] = dict(self.present_env)
            outputdict['elements_low_entropy']=self.list_to_remove
            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)


            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling


            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                       0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                       0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.Plot_PSE_DICT = {}
        self.present_env = {}
        # valence dependency can be introduced later
        # counter_not_primitive=0
        for mat in list_mat:
            #print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)
            # prim_struct=lse.structure.get_primitive_structure()
            # print(prim_struct.composition)
            # print(lse.structure.composition)
            # if prim_struct.composition==lse.structure.composition:


            pauling0 = Pauling0(lse)
            pauling5 = Pauling5()
            if not self.start_from_connections or not os.path.isfile(
                    os.path.join(self.connections_folder, mat + '.json')):
                pauling5.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder, distance=8.0)
            else:
                pauling5.from_file(filename=mat + '.json', foldername=self.connections_folder)

            try:
                if pauling5.is_fulfilled(leave_out_list=self.list_to_remove):
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)

            try:
                Details = pauling5.get_details(leave_out_list=self.list_to_remove)
                #print(Details)
                New_Details = self._reformat_details_elementdependency(Details)

                self._add_dict_cat_dependency(self.Plot_PSE_DICT, New_Details)
                self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                              number_of_elements_to_add=1)

            except RuleCannotBeAnalyzedError:
                pass

            # else:
            #     counter_not_primitive+=1
            #     print(counter_not_primitive)

    def _fifth_rule_plot(self, okay, notokay, title="Tested Structures"):
        import matplotlib.pyplot as plt

        labels = 'Fulfills rule', 'Does not fulfill rule'
        sizes = [okay, notokay]
        colors = ['#5da5daff', '#4d4d4dff']
        explode = (0, 0)  # explode 1st slice

        # Plot
        plt.pie(sizes, explode=explode, labels=labels, colors=colors,
                autopct='%1.1f%%', shadow=False, startangle=140)
        plt.title(title)
        plt.axis('equal')
        return plt

    def _reformat_details_elementdependency(self, details):
        new_dict = {}
        for key in details:
            if not key in new_dict:
                new_dict[key] = [details[key]['fulfilled'], details[key]['not_fulfilled']]
        return new_dict


class AllPaulingOverAllAnalysis(OverAllAnalysis):
    def run(self, remove_elements_low_entropy=False, start_from_connections=False,
            save_connections=True, connections_folder34='AnalysisConnections', connections_folder5='AnalysisConnections_5thRule',
            start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analyisis=False, save_structure_analysis=True,
            path_to_save='Results/Results_AllRules.json',threshold_remove_elements=0.95,start_material=None,stop_material=None):
        self.start_material = start_material
        self.stop_material = stop_material

        self.remove_elements_low_entropy = remove_elements_low_entropy
        self.connections_folder34 = connections_folder34
        self.connections_folder5 = connections_folder5
        self.start_from_connections = start_from_connections
        self.save_connections = save_connections
        self.threshold_remove_elements=threshold_remove_elements


        if not start_from_results:
            if self.remove_elements_low_entropy:
                     newclass = Pauling1Entropy(source=self.source, onlybinaries=self.onlybinaries,
                                                plot_element_dependend_analysis=False,list_of_materials_to_investigate=self.list_of_materials_to_investigate,start_material=self.start_material,stop_material=self.stop_material)
                     newclass.run(start_from_results=False, save_result_data=False)
                     self.list_to_remove = [el for el, item in newclass.Plot_PSE_entropy.items() if
                                            item > self.threshold_remove_elements]
            else:
                self.list_to_remove=[]

            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.list_to_remove=inputdict['list_to_remove']

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['list_to_remove']=self.list_to_remove
            self._save_results_to_file(outputdict,path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analyisis,
                                                                             start_from_Matching=self.use_prematching)

        if save_structure_analysis:
            self._print_to_file_similar_structures(dict_similarstructures_exceptions, filename=path_to_save.split('.')[
                                                                                                   0] + "_structural_exceptions_readable.yaml")

            self._print_to_file_similar_structures(dict_similarstructures_fulfilling, filename=path_to_save.split('.')[
                                                                                                   0] + "_structures_fulfilling_readable.yaml")

            self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                   filename=path_to_save.split('.')[
                                                                0] + "_structural_exceptions_readable.csv")

            self._print_to_file_similar_structures(dict_similarstructures_fulfilling, fmt='csv',
                                                   filename=path_to_save.split('.')[
                                                                0] + "_structures_fulfilling_readable.csv")

    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        # how to do an elementwise analysis?
        # self.Dict_val = {}
        # self.Dict_CN = {}
        # self.Plot_PSE_DICT = {}
        # self.present_env = {}
        # valence dependency can be introduced later
        #print(self.list_to_remove)
        for mat in list_mat:
            lse = self._get_lse_from_folder(mat, source=self.source)
            pauling1 = Pauling1(lse, filenameradii='../univalent_cat_radii.json', onlylowerlimit=False)
            pauling2 = Pauling2(lse)
            pauling3 = Pauling3()
            pauling4 = Pauling4()
            pauling5 = Pauling5()
            if not self.start_from_connections or not os.path.isfile(
                    os.path.join(self.connections_folder34, mat + '.json')) or not os.path.isfile(
                os.path.join(self.connections_folder5, mat + '.json')):
                pauling3.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder34, distance=8.0)
                pauling4.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder34, distance=8.0)
                pauling5.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder5, distance=8.0)
            else:
                pauling3.from_file(filename=mat + '.json', foldername=self.connections_folder34)
                pauling4.from_file(filename=mat + '.json', foldername=self.connections_folder34)
                pauling5.from_file(filename=mat + '.json', foldername=self.connections_folder5)

            # pauling0 = Pauling0(lse)
            try:
                result1=pauling1.is_fulfilled()
                result2=pauling2.is_fulfilled()
                result3=pauling3.is_fulfilled()
                result4=pauling4.is_fulfilled()
                result5=pauling5.is_fulfilled(leave_out_list=self.list_to_remove)

                if result1 and result2 and result3 and result4 and result5:
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)

            # How to analyse elements in more detail?s

            #
            # try:
            #     Details = pauling5.get_details(leave_out_list=self.list_to_remove)
            #     New_Details = self._reformat_details_elementdependency(Details)
            #     # print(New_Details)
            #
            #     self._add_dict_cat_dependency(self.Plot_PSE_DICT, New_Details)
            #     # !!!!! not working correcntly
            #     print("Problem")
            #     self._add_dict_cat_dependency(self.present_env, pauling0._get_cations_in_structure(),
            #                                   number_of_elements_to_add=1)
            # #
            # except RuleCannotBeAnalyzedError:
            #     print("Exception")

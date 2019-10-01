import json
import numpy as np
from collections import OrderedDict, Counter
import matplotlib
import matplotlib.pyplot as plt
import os
from yaml import dump
import csv

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.chemenv.utils.chemenv_errors import NeighborsNotComputedChemenvError
from pymatgen.analysis.chemenv.utils.scripts_utils import draw_cg
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import AllCoordinationGeometries
from pymatgen.core.sites import PeriodicSite
from pymatgen.vis.structure_vtk import StructureVis
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
from pymatgen.core.periodic_table import Element

from PaulingRules import Pauling0, Pauling1, Pauling2, Pauling3, Pauling4, Pauling5, RuleCannotBeAnalyzedError, \
    FrequencyEnvironmentPauling1, get_entropy_from_frequencies, get_most_frequent_environment, \
    get_mean_CN_from_frequencies
from PlotClasses import PlotterPSE


class OverAllAnalysis:

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 analyse_structures=True, use_prematching=True, list_of_materials_to_investigate=None,print_structure_matching=True):
        """

        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :param onlybinaries: only binary structures will be analysed
        :param plot_element_dependend_analysis: will plot element dependent anlysis
        :param lowest_number_environments_for_plot: decides which elements will be considered in element dependent plot with the help of the lowest number of environments
        :param lower_limit_plot: decides the lower limit for plotting this
        :param upper_limit_plot: decides the upper limit for plotting this
        :param analyse_structures: should the structures be analysed and matched
        :param use_prematching: decides whether a already existing matching will be used, only available for data from Materials Project
        :param list_of_materials_to_investigate: you can also a reference to a list of materials that should be investigated, e.g. "test.json"
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
        self.print_structure_matching=print_structure_matching
        # could include a function that starts plotting from saved data?
        # how should one do this in the best way

    def _get_list_materials(self, source='MP', onlybinaries=False, start_material=None, stop_material=None) -> list:
        """
        will read lists of materials from file
        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :param onlybinaries: from the list of materials, only binaries are used
        :param start_material: number that cuts the list of all materials
        :param stop_material: number that cuts the list of all materials
        :return: list of materials
        """
        if source == 'MP':
            with open("../Assessment/Should_not_be_changed/allmaterials.json", "r") as f:
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
                    list_compound = list_compound_dict["is_clear_compounds"]

        elif source == 'MP_very_symmetric':
            with open(
                    "../Assessment/Should_not_be_changed/ce_fraction_0.95plus_csm_0.1plus_eabovehull0.025plus_discardS1.json",
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
                    list_compound = list_compound_dict["is_clear_compounds"]

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

        elif source == 'experimental':
            # TODO: update
            list_compound = self._get_precomputed_results(
                "../Assessment/Should_not_be_changed/List_experimental_oxides.json")
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

    def _get_lse_from_folder(self, mat: str, source='MP') -> LightStructureEnvironments:
        """
        will fetch lse from a folder
        :param mat: name of the material, e.g. "mp-7000"
        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :return: LightStructureEnvironments
        """
        if source == 'MP' or source == 'MP_very_symmetric' or source == 'my_own_list':
            with open(os.path.join("../Assessment/lse_MP", mat + ".json"), 'r') as f:
                data = json.load(f)
        elif source == 'experimental':
            with open(os.path.join("../Assessment/lse_exp", mat + ".json"), 'r') as f:
                data = json.load(f)

        lse = LightStructureEnvironments.from_dict(data)
        return lse

    def _plot_PSE(self, Dict_to_Plot: dict, lowest_number_of_environments_considered: int,
                  upper_number_of_environments_considered=None, xlim=[1, 18], ylim=[1, 10],
                  lowerlimit=0,
                  upperlimit=1, counter_cations_env=None, plot_directly_from_freq=False) -> plt:
        """
        will plot data in a periodic table
        :param Dict_to_Plot: dict that will be plotted
        :param lowest_number_of_environments_considered:
        :param upper_number_of_environments_considered:
        :param xlim: xrange for periodic table
        :param ylim: yrange for periodic table
        :param lowerlimit: lower limit for plot
        :param upperlimit: upper limit for plot
        :param counter_cations_env: dict that tells you how many environments are present for each cation
        :param plot_directly_from_freq: will directly plot the entry for each element in the dict
        :return: plot
        """

        plotterpse = PlotterPSE(valuestoplot=Dict_to_Plot, counter_cations_env=counter_cations_env,
                                plot_directly_from_freq=plot_directly_from_freq)
        # atoms to consider for plot instead?

        plt = plotterpse.get_plot(xlim=xlim, ylim=ylim,
                                  lowest_number_of_environments_considered=lowest_number_of_environments_considered,
                                  upper_number_of_environments_considered=upper_number_of_environments_considered,
                                  lowerlimit=lowerlimit, upperlimit=upperlimit)
        return plt

    def _save_results_to_file(self, dict_to_save, path):
        with open(path, 'w') as f:
            json.dump(dict_to_save, f)

    def _get_precomputed_results(self, path):
        with open(path, 'r') as f:
            dict_here = json.load(f)
        return dict_here

    def _add_dict_cat_dependency(self, start_dict: dict, dict_to_add: dict, number_of_elements_to_add=2):
        """

        :param start_dict: dict that will be extended
        :param dict_to_add: dict that will be used to extend
        :param number_of_elements_to_add:
             ==1: following dicts can be summed, integers for each key will be summed: {"Ga":1, "Sn":2}
             ==2: following dicts can be summed, individual numbers in the arrays of each key will be summed: {"Ga": [1,0], "Sn": [0,2]}
             ==4: following dicts can be extended with another dict: {"Ga":["O:6","O:6"], "Sn":["T:4","T:4"]}
        :return: None, start_dict will have the new values
        """

        for key, item in dict_to_add.items():
            if not key in start_dict:
                if number_of_elements_to_add == 1:
                    start_dict[key] = dict_to_add[key]
                elif number_of_elements_to_add == 2:
                    start_dict[key] = dict_to_add[key].copy()
                elif number_of_elements_to_add == 4:
                    start_dict[key] = dict_to_add[key].copy()
            else:
                if number_of_elements_to_add == 1:
                    start_dict[key] += dict_to_add[key]
                elif number_of_elements_to_add == 2:
                    start_dict[key][0] += dict_to_add[key][0]
                    start_dict[key][1] += dict_to_add[key][1]
                elif number_of_elements_to_add == 4:
                    start_dict[key].extend(dict_to_add[key])

    def _get_similar_structures(self, list_mat_id: list, source='MP', save_to_file=True,
                                path_to_save='Similar_Structures.json', fetch_results_only=False,
                                start_from_Matching=False,
                                path_to_precomputed_matching="Should_not_be_changed/Matching_All_Structures.json",
                                restart_from_matching=False) -> dict:
        """
        will match materials ids according to structures
        :param list_mat_id: list of materials ids
        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :param save_to_file:
        :param path_to_save:
        :param fetch_results_only:
        :param start_from_Matching: matches with the help of existing matching
        :param restart_from_matching: continues a matching from a file, e.g. to continue it later when interrupted
        :return: dict including the matching
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
                    # TODO: check if this is okay! maybe a different algorithm is needed to do so?
                    # TODO: maybe test function from pymatgen
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
            if source == 'MP' or source == 'MP_very_symmetric' or source == 'my_own_list':
                prematching = self._get_precomputed_results(path_to_precomputed_matching)
            elif source == 'experimental':
                Warning.warns("No pre-matching exists")
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
                        # TODO: check if this can ever happen
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
        # TODO: test valueerror -> use another list before loading (start and stop_material nutzen)
        if fetch_results_only:
            outputdict = self._get_precomputed_results(path_to_save)
            if not set(outputdict['list_mat_id']) == set(list_mat_id):
                raise ValueError
        return outputdict

    def _print_to_file_similar_structures(self, dict_to_print: dict, source='MP', filename='Test.yaml', fmt='yml',
                                          add_info=None, name_add_info=""):
        """
        will print the matched structures in a file
        :param dict_to_print: dict to print
        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :param filename: filename
        :param fmt: 'yml' or 'csv'
        :param add_info: additional external output for files
        :param name_add_info: how is the external output called

        """

        if self.print_structure_matching:
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

                        if add_info is None:
                            new_dict_to_print[key][item] = {"formula": dict_to_print['additional_info'][item], "CN": CN,
                                                            "valences": valence, "cations": cat}
                        else:
                            new_dict_to_print[key][item] = {"formula": dict_to_print['additional_info'][item], "CN": CN,
                                                            "valences": valence, "cations": cat,
                                                            name_add_info: add_info[item]}

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
                        if add_info is None:
                            new_dict_to_print.append(
                                {"mpid": str(item), "formula": dict_to_print['additional_info'][item], "CN": CN,
                                 "valences": valence, "cations": cat, "structure_type": key})
                        else:
                            new_dict_to_print.append(
                                {"mpid": str(item), "formula": dict_to_print['additional_info'][item], "CN": CN,
                                 "valences": valence, "cations": cat, "structure_type": key, name_add_info: add_info[item]})

                if add_info is None:
                    with open(filename, 'w') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter=';',
                                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                        filewriter.writerow(['mp-id', 'formula', 'structure_type', 'cations', 'valences', 'CNs'])
                        for line in new_dict_to_print:
                            filewriter.writerow(
                                [str(line["mpid"]), str(line["formula"]), str(line["structure_type"]), str(line["cations"]),
                                 str(line["valences"]), str(line["CN"])])
                else:
                    with open(filename, 'w') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter=';',
                                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                        filewriter.writerow(
                            ['mp-id', 'formula', 'structure_type', 'cations', 'valences', 'CNs', name_add_info])
                        for line in new_dict_to_print:
                            filewriter.writerow(
                                [str(line["mpid"]), str(line["formula"]), str(line["structure_type"]), str(line["cations"]),
                                 str(line["valences"]), str(line["CN"]), str(line[name_add_info])])

        else:
            if fmt == 'yml':
                new_dict_to_print = OrderedDict()
                for key, items in OrderedDict(
                        sorted(dict_to_print['structure_matching'].items(), key=lambda t: len(t[1]),
                               reverse=True)).items():

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

                        if add_info is None:
                            new_dict_to_print[key][item] = {"formula": dict_to_print['additional_info'][item], "CN": CN,
                                                            "valences": valence, "cations": cat}
                        else:
                            new_dict_to_print[key][item] = {"formula": dict_to_print['additional_info'][item], "CN": CN,
                                                            "valences": valence, "cations": cat,
                                                            name_add_info: add_info[item]}

                with open(filename, 'w') as f:
                    dump(OrderedDict(new_dict_to_print), f)

            elif fmt == 'csv':
                new_dict_to_print = []
                for key, items in OrderedDict(
                        sorted(dict_to_print['structure_matching'].items(), key=lambda t: len(t[1]),
                               reverse=True)).items():

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
                        if add_info is None:
                            new_dict_to_print.append(
                                {"mpid": str(item), "formula": dict_to_print['additional_info'][item], "CN": CN,
                                 "valences": valence, "cations": cat})
                        else:
                            new_dict_to_print.append(
                                {"mpid": str(item), "formula": dict_to_print['additional_info'][item], "CN": CN,
                                 "valences": valence, "cations": cat,
                                 name_add_info: add_info[item]})

                if add_info is None:
                    with open(filename, 'w') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter=';',
                                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                        filewriter.writerow(['mp-id', 'formula', 'cations', 'valences', 'CNs'])
                        for line in new_dict_to_print:
                            filewriter.writerow(
                                [str(line["mpid"]), str(line["formula"]),
                                 str(line["cations"]),
                                 str(line["valences"]), str(line["CN"])])
                else:
                    with open(filename, 'w') as csvfile:
                        filewriter = csv.writer(csvfile, delimiter=';',
                                                quotechar='|', quoting=csv.QUOTE_MINIMAL)
                        filewriter.writerow(
                            ['mp-id', 'formula', 'cations', 'valences', 'CNs', name_add_info])
                        for line in new_dict_to_print:
                            filewriter.writerow(
                                [str(line["mpid"]), str(line["formula"]),
                                 str(line["cations"]),
                                 str(line["valences"]), str(line["CN"]), str(line[name_add_info])])

    def _pieplot_connections(self, corner: int, edge: int, face: int, title="All") -> plt:
        """
        will plot a pie chart based on integers describing the connections
        :param corner: number of corner connections
        :param edge: number of edge connections
        :param face: number of face connections
        :param title: title of the plot
        :return: plot
        """
        labels = 'via corners', 'via edges', 'via faces'
        sizes = [corner, edge, face]
        colors = ['#5da5daff', '#faa43aff', '#f15854ff']

        plt.pie(sizes, labels=labels, colors=colors,
                autopct='%1.1f%%', shadow=False, startangle=140)
        plt.title(title)
        plt.axis('equal')
        return plt

    def visualize_structure_by_id(self, mat, supercell=[1, 1, 0]):
        """

        :param mat: materials id
        :param supercell: list of additional cells that are considered in each direction
        :return:
        """
        # TODO: rethink supercell implementation
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


class HowMany(OverAllAnalysis):
    """Class to analyse the largest frequency of a coordination environment for each cation"""

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lower_limit_plot=0.0, upper_limit_plot=1.0,
                 list_of_materials_to_investigate=None, start_material=None, stop_material=None,
                 lowest_number_of_environments_considered=0, upper_number_of_environments_considered=None):
        """
            :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
            :param onlybinaries: only binary structures will be analysed
            :param save_result_data: all results will be saved so that results can easily be recreated
            :return:
        """

        self.source = source
        self.onlybinaries = onlybinaries
        self.plot_element_dependend_analysis = plot_element_dependend_analysis
        self.lower_limit_plot = lower_limit_plot
        self.upper_limit_plot = upper_limit_plot
        self.lowest_number_of_environments_considered = lowest_number_of_environments_considered
        self.upper_number_of_environments_considered = upper_number_of_environments_considered
        self.start_material = start_material
        self.stop_material = stop_material
        self.list_of_materials_to_investigate = list_of_materials_to_investigate

    def run(self, start_from_results=False, save_result_data=True,
            path_to_save='Results/Results_Number_Compounds.json'):
        """
        runs the analysis
        :param start_from_results: if True, restart from Result file
        :param save_result_data: if True, saves results in json file
        :param path_to_save: path to save the files
        :return:
        """
        if not start_from_results:

            self._new_setup()

        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.present_env = Counter(inputdict['Counter_cation'])

        if self.plot_element_dependend_analysis:
            # TODO put a parameter that can exclude certain elements!
            # TODO: test this class
            plt = self._plot_PSE(self.present_env,
                                 plot_directly_from_freq=True, lowerlimit=self.lower_limit_plot,
                                 lowest_number_of_environments_considered=self.lowest_number_of_environments_considered,
                                 upperlimit=self.upper_limit_plot,
                                 upper_number_of_environments_considered=self.upper_number_of_environments_considered,
                                 counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['Counter_cation'] = dict(self.present_env)
            self._save_results_to_file(outputdict, path_to_save)

    def _new_setup(self):
        """
        could count the number of environments instead of structures?
        will be done if results are not read from file
        :return:
        """
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)
        self.All_Details = {}
        self.present_env = {}
        for mat in list_mat:
            # print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)
            pauling0 = Pauling0(lse)
            self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                          number_of_elements_to_add=1)


class Pauling1Frequency(OverAllAnalysis):
    """Class to analyse the largest frequency of a coordination environment for each cation"""

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 list_of_materials_to_investigate=None, start_material=None, stop_material=None):
        """
            :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
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
        self.start_material = start_material
        self.stop_material = stop_material
        self.list_of_materials_to_investigate = list_of_materials_to_investigate

    def run(self, start_from_results=False, save_result_data=True, path_to_save='Results/Results_First_Limits.json'):
        """
        runs the analysis
        :param start_from_results: if True, restart from Result file
        :param save_result_data: if True, saves results in json file
        :param path_to_save: path to save the files
        :return:
        """
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
        """
        will be done if results are not read from file
        :return:
        """
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)
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
    """Class to analyse the Shannon entropy for each element"""

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 list_of_materials_to_investigate=None, start_material=None, stop_material=None):
        """

        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :param onlybinaries: only binary structures will be analysed
        :param plot_element_dependend_analysis: will plot element dependent anlysis
        :param lowest_number_environments_for_plot: decides which elements will be considered in element dependent plot with the help of the lowest number of environments
        :param lower_limit_plot: decides the lower limit for plotting this
        :param upper_limit_plot: decides the upper limit for plotting this
        :param list_of_materials_to_investigate: you can also a reference to a list of materials that should be investigated, e.g. "test.json"
        :param start_material: number at which the investigation is started
        :param stop_material: number before which the investigation is started
        """

        super(Pauling1Entropy, self).__init__(source=source, onlybinaries=onlybinaries,
                                              plot_element_dependend_analysis=plot_element_dependend_analysis,
                                              lowest_number_environments_for_plot=lowest_number_environments_for_plot,
                                              lower_limit_plot=lower_limit_plot, upper_limit_plot=upper_limit_plot,
                                              list_of_materials_to_investigate=list_of_materials_to_investigate,
                                              start_material=start_material, stop_material=stop_material)

    def run(self, start_from_results=False, save_result_data=True, path_to_save='Results/Results_First_Limits.json'):
        """

        :param start_from_results: if True, restart from Result file
        :param save_result_data: if True, saves results in json file
        :param path_to_save: path to save the files
        :return:
        """
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


class Pauling1MeanCoordinationNumber(Pauling1Entropy):
    """
    Class to analyse the mean CN for each element
    """

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=2.0, upper_limit_plot=12,
                 list_of_materials_to_investigate=None, start_material=None, stop_material=None):
        """

        :param source: 'MP' (Materials Project), 'MP_very_symmetric' (only structures with very symmetric coordiation environments), or 'experimental' (structures from COD) can be used
        :param onlybinaries: only binary structures will be analysed
        :param plot_element_dependend_analysis: will plot element dependent anlysis
        :param lowest_number_environments_for_plot: decides which elements will be considered in element dependent plot with the help of the lowest number of environments
        :param lower_limit_plot: decides the lower limit for plotting this
        :param upper_limit_plot: decides the upper limit for plotting this
        :param list_of_materials_to_investigate: you can also a reference to a list of materials that should be investigated, e.g. "test.json"
        :param start_material: number at which the investigation is started
        :param stop_material: number before which the investigation is started
        """
        super(Pauling1Entropy, self).__init__(source=source, onlybinaries=onlybinaries,
                                              plot_element_dependend_analysis=plot_element_dependend_analysis,
                                              lowest_number_environments_for_plot=lowest_number_environments_for_plot,
                                              lower_limit_plot=lower_limit_plot, upper_limit_plot=upper_limit_plot,
                                              list_of_materials_to_investigate=list_of_materials_to_investigate,
                                              start_material=start_material, stop_material=stop_material)

    def run(self, start_from_results=False, save_result_data=True, path_to_save='Results/Results_First_Limits.json'):
        """
        will run the analysis
        :param start_from_results: if True, restart from Result file
        :param save_result_data: if True, saves results in json file
        :param path_to_save: path to save the files
        :return:
        """

        if not start_from_results:
            super(Pauling1Entropy, self)._new_setup()
            self.Plot_PSE_numbers = get_mean_CN_from_frequencies(self.All_Details)
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.Plot_PSE_numbersE = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])

        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_numbers,
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 plot_directly_from_freq=True, lowerlimit=self.lower_limit_plot,
                                 upperlimit=self.upper_limit_plot, counter_cations_env=self.present_env)
            plt.show()

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['PSE_Dict'] = self.Plot_PSE_numbers
            outputdict['Counter_cation'] = dict(self.present_env)
            self._save_results_to_file(outputdict, path_to_save)


class Pauling1OverAllAnalysis(OverAllAnalysis):
    """
    Class to analyse the first rule based on univalent radii
    """

    def run(self, start_from_results=False, save_result_data=True, save_structure_analysis=True,
            restart_from_saved_structure_analysis=False,
            path_to_save='Results/Results_First_Rule.json', start_material=None, stop_material=None):
        """

        :param start_from_results: will start from file if true
        :param save_result_data: will save a file
        :param save_structure_analysis: will save the structure analysis as well (path will be set automatically)
        :param restart_from_saved_structure_analysis: will start from a save structure analysis (to continue if stops)
        :param path_to_save: path at which the normal results are saved
        :param start_material: value cuts the list of all materials
        :param stop_material: value cuts the list of all materials
        :return:
        """
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

        # print result to screen:
        all = 0
        fulfilling = 0
        for value in self.Plot_PSE_DICT.values():
            fulfilling += value[0]
            all += value[0]
            all += value[1]
        print("Only " + str(float(fulfilling) / float(all)) + " of all environments fulfill the first rule.")
        print("This many environments were considered in the analysis:" + str(all))
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
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)
            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split(
                                                                                 '.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling

            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _new_setup(self):
        """
        will do the analysis from scratch
        :return:
        """
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

            pauling1 = Pauling1(lse=lse, filenameradii='../Assessment/Should_not_be_changed/univalent_cat_radii.json',
                                onlylowerlimit=False)
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


class EntropyDeviationFrom2ndRuleDiagram(OverAllAnalysis):
    """
    Class to find explanation for bad performance of the second rule

    """

    def run(self, show_plot=True, start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analysis=False,
            save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule_Entropy_vs_Deviation.json', start_material=None,
            stop_material=None):
        # TODO: will start first second rule and combine both to make a plot
        # several options should be considered: avg, min, max entropy and av, min, max deviation
        #get additional_info from 2nd rule
        #get entropies from first rule

        #then loop over all  materials get avg, min, max entropy, and avg, min, max deviation

        #make 9 plots and check if there is any correlation
        self.start_material=start_material
        self.stop_material=stop_material

        newclass = Pauling1Entropy(source=self.source, onlybinaries=self.onlybinaries,
                                   plot_element_dependend_analysis=False,
                                   list_of_materials_to_investigate=None,
                                   start_material=self.start_material, stop_material=self.stop_material)
        newclass.run(start_from_results=False, save_result_data=False)
        entropy=newclass.Plot_PSE_entropy

        newclass2 = Pauling2OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=False,
                                           lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=0.8,
                                           analyse_structures=False, use_prematching=True)
        newclass2.run(start_from_results=False, save_result_data=False, path_to_save='Results/Results_Second_Rule.json',
                     save_structure_analysis=True, restart_from_saved_structure_analysis=False, show_histogram=False, stepsize_histogram=0.1,show_plot=False)

        self.additional_info=newclass2.additional_info

        #brauche composition der struktur
        #daraus dann ermitteln, was die entropie ist (avg, min,max)
        #gegen alle moeglichen deviations (avg,min,max) plotten


        #TODO: laufe ueber alle mps und strukturen und sammle alle
        min_dev=[]
        max_dev=[]
        mean_dev=[] #think of another way, will always be zero
        min_entropy_list=[]
        max_entropy_list=[]
        avg_entropy_list=[]
        for key,value in self.additional_info.items():
            min_dev.append(min([abs(x-2.0) for x in value]))
            max_dev.append(max([abs(x-2.0) for x in value]))
            mean_dev.append(np.mean([abs(x-2.0) for x in value]))

            lse = self._get_lse_from_folder(mat=key, source=self.source)
            composition=lse.structure.composition
            # print(composition)

            elements=composition.elements
            minentropy=1.0
            maxentropy=0.0
            avgenentropy_sum=0.0
            nb_cations=0
            for el in elements:
                #print(el)
                if str(el)!='O':
                    test=entropy[str(el)]
                    nb_cations+=1
                    avgenentropy_sum+=test
                    #print(test)
                    if test <= minentropy:
                        minentropy=test
                    if test >= maxentropy:
                        maxentropy=test
                    #print(minentropy)
                    #print(maxentropy)
            #exit()

            max_entropy_list.append(maxentropy)
            min_entropy_list.append(minentropy)
            avg_entropy_list.append(avgenentropy_sum/nb_cations)

        plt.plot(max_entropy_list,min_dev,'x')
        plt.show()
        plt.plot(max_entropy_list,max_dev,'x')
        plt.show()
        plt.plot(max_entropy_list,mean_dev,'x')
        plt.show()


        plt.plot(min_entropy_list,min_dev,'x')
        plt.show()
        plt.plot(min_entropy_list,max_dev,'x')
        plt.show()
        plt.plot(min_entropy_list,mean_dev,'x')
        plt.show()

        plt.plot(avg_entropy_list,min_dev,'x')
        plt.show()
        plt.plot(avg_entropy_list,max_dev,'x')
        plt.show()
        plt.plot(avg_entropy_list,mean_dev,'x')
        plt.show()



class Pauling2OverAllAnalysis(OverAllAnalysis):
    """
    Class to analyse second rule
    """

    def run(self, show_plot=True, start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analysis=False,
            save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json', start_material=None, stop_material=None,
            show_histogram=False, stepsize_histogram=0.1):
        """
        :param show_plot: will show the main analysis plot from the second rule
        :param start_from_results: if True,   restart from saved results
        :param save_result_data: if True, save results in file
        :param restart_from_saved_structure_analysis: if True, restarts a started structure analysis
        :param save_structure_analysis:  If True: saves structure analysis
        :param path_to_save: path to save normal results, structure analysis will saved in an adapted path
        :param start_material: number of material at which the analysis starts
        :param stop_material: number of material before which the analysis stops
        :param show_histogram: if True, shows histogram with +- deviations for second rule
        :param stepsize_histogram: how large are the charge steps for the histogram
        :return:
        """

        if not start_from_results:
            self.start_material = start_material
            self.stop_material = stop_material
            self.stepsize_histogram = stepsize_histogram
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
            self.additional_info = inputdict["additional_info"]
            self.arraydev_share_plus_minus = np.array(inputdict['arraydev_share_plus_minus'])
            self.relativefrequency_plus_minus = np.array(inputdict['relativefrequency_plus_minus'])
            self.stepsize_histogram = inputdict['stepsize_histogram']
            self.extreme_exceptions = inputdict['list_extreme_exceptions']

        if show_plot:
            plot = self._secondrule_plot(arraydev=self.arraydev_share * 100.0,
                                         relativefreqarray=self.relativefrequency * 100.0,
                                         tot_stddev=self.tot_stddev * 100.0)
            plot.show()

        if show_histogram:
            plt2 = self._secondrule_plot_histogram(self.arraydev_share_plus_minus * 100.0,
                                                   self.relativefrequency_plus_minus * 100.0)
            plt2.show()

        if self.plot_element_dependend_analysis:
            plt = self._plot_PSE(self.Plot_PSE_DICT, xlim=[1, 18], ylim=[1, 10],
                                 lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                 lowerlimit=self.lower_limit_plot, upperlimit=self.upper_limit_plot,
                                 counter_cations_env=self.present_env)
            plt.show()

        # TODO: test this part!
        self.extreme_exceptions = self._get_extreme_exceptions(self.additional_info, perc_std=self.tot_stddev[0],
                                                               ideal_value=2, larger_how_many_std=3)

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
            outputdict['additional_info'] = self.additional_info
            outputdict['arraydev_share_plus_minus'] = list(self.arraydev_share_plus_minus)
            outputdict['relativefrequency_plus_minus'] = list(self.relativefrequency_plus_minus)
            outputdict['stepsize_histogram'] = self.stepsize_histogram
            outputdict['list_extreme_exceptions'] = self.extreme_exceptions

            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             start_from_Matching=self.use_prematching,
                                                                             fetch_results_only=restart_from_saved_structure_analysis)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             start_from_Matching=self.use_prematching,
                                                                             fetch_results_only=restart_from_saved_structure_analysis)

            dict_similarstructures_extreme_exceptions = self._get_similar_structures(self.extreme_exceptions,
                                                                                     source=self.source,
                                                                                     save_to_file=save_structure_analysis,
                                                                                     path_to_save=
                                                                                     path_to_save.split('.')[
                                                                                         0] + "_structural_extreme_exceptions.json",
                                                                                     start_from_Matching=self.use_prematching,
                                                                                     fetch_results_only=restart_from_saved_structure_analysis)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling
            self.dict_similarstructures_extreme_exceptions = dict_similarstructures_extreme_exceptions

            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml",
                                                       add_info=self.additional_info, name_add_info="Bond valence sum")

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml",
                                                       add_info=self.additional_info, name_add_info="Bond valence sum")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv",
                                                       add_info=self.additional_info, name_add_info="Bond valence sum")

                self._print_to_file_similar_structures(dict_similarstructures_extreme_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_extreme_exceptions_readable.csv",
                                                       add_info=self.additional_info, name_add_info="Bond valence sum")

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv",
                                                       add_info=self.additional_info, name_add_info="Bond valence sum")

    def _get_extreme_exceptions(self, additional_info: dict, perc_std: float, ideal_value: float,
                                larger_how_many_std=2.0) -> list:
        """
        filters for materials names with extreme exceptions
        :param additional_info: dict {"materialname":[2.3,2,0.8]
        :param perc_std: value between 0 and 1
        :param ideal_value: usually 2 for oxides
        :param larger_how_many_std: how many std deviations should be considered
        :return:
        """

        new_list = []
        for key, value in additional_info.items():
            if (max([abs(float(ideal_value) - float(x)) for x in value]) - float(larger_how_many_std) * float(
                    perc_std) * float(ideal_value)) > 10e-8:
                new_list.append(key)
        return new_list

    def _stddev(self, lst: list, mean: float) -> list:
        """
        calculates standard deviation in percent of a sample in a list format

        :param lst: lst for which the standard deviation is calculated
        :param mean: float describing the mean of the lst
        :return: list of element including the standard deviation
        """

        # TODO: get rid of the list format
        sum = 0.0
        mn = mean
        for i in range(len(lst)):
            sum += pow((lst[i] - mn), 2)
        return np.sqrt([sum / (len(lst) - 1)]) / mn

    def _get_absolute_deviation_from_ideal_value(self, inputarray: list, ideal: float) -> list:
        """
        calculates the deviation from an ideal value for a whole list
        :param inputarray:
        :param ideal:
        :return:
        """
        return [abs(float(x) - float(ideal)) for x in inputarray]

    def _get_deviation_from_ideal_value(self, inputarray: list, ideal: float) -> list:
        """
        calculates the deviation from an ideal value for a whole list
        :param inputarray:
        :param ideal:
        :return:
        """
        # TODO: test
        return [(float(x) - float(ideal)) for x in inputarray]

    def _get_frequency_of_values(self, arraydeviations: list, stepsize: float) -> tuple:
        """
        will get the frequencies of values lower a certain step size
        :param arraydeviations: array of the deviations from ideal value
        :param stepsize: stepsize as a float
        :return:
        """
        frequency = []
        dev_array = np.arange(0, max(arraydeviations) + stepsize, stepsize)

        for step in dev_array:
            # print(step)

            frequency.append(len([x for x in arraydeviations if x <= step]))
            # print(frequency)

        return dev_array, frequency

    def _get_frequency_bs_for_plus_minus_deviation(self, arraydeviations: list, stepsize: float) -> tuple:
        """
        will get the frequencies of values lower a certain step size
        :param arraydeviations: array of the deviations from ideal value
        :param stepsize: stepsize as a float
        :return:
        """
        frequency = []
        dev_array = np.arange(min(arraydeviations), max(arraydeviations) + stepsize, stepsize)

        for step in dev_array:
            frequency.append(
                len([x for x in arraydeviations if x <= step + 0.5 * stepsize and x > (step - 0.5 * stepsize)]))

        # print(dev_array)
        # print(frequency)
        # print(sum(frequency))
        # print(len(arraydeviations))
        return dev_array, frequency

    def _list_to_np_array_and_divide_by_value(self, inputarray: list, valuetodivide: float) -> np.array:
        """
        will divide each element of a whole array by a value
        :param inputarray: array
        :param valuetodivide: float
        :return: numpy.array
        """
        return (np.array(inputarray) / float(valuetodivide))

    def _secondrule_plot_histogram(self, arraydev, relativefreq):
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        font = {'size': 22}

        matplotlib.rc('font', **font)
        plt.bar(self.arraydev_share_plus_minus * 100, height=self.relativefrequency_plus_minus * 100)
        plt.xlabel("Deviation (%) from the ideal valence -2")
        plt.ylabel("Oxygen Atoms (%)")

        return plt

    def _secondrule_plot(self, arraydev: np.array, relativefreqarray: np.array, tot_stddev: float, maxpercentage=70,
                         save_plot=False,
                         filename='second_rule.svg') -> plt:
        """

        :param arraydev: will be plotted in x direction
        :param relativefreqarray: will be plotted in y direction
        :param tot_stddev: total standard deviation in percent
        :param maxpercentage: up to which percentage should the plot be shown
        :param save_plot: should the plot be saved
        :param filename: filename for saving this plot
        :return:
        """
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
        """
        does the analysis from scratch
        :return:
        """
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)
        # print("Number of Materials")
        # print(len(list_mat))
        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        array_bvs = []
        self.Plot_PSE_DICT = {}
        self.present_env = {}
        self.additional_info = {}
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
            # save a dict that can be used later for printing additional info
            # TODO: has to be tested
            self.additional_info[mat] = bvs

        self.array_bvs = array_bvs
        self.bs_sum_mean = np.mean(array_bvs)

        self.tot_stddev = self._stddev(array_bvs, np.mean(array_bvs))
        ideal_bs = 2.0
        self.ideal_bs = 2.0
        bs_dev = self._get_absolute_deviation_from_ideal_value(inputarray=array_bvs, ideal=ideal_bs)
        self.bs_dev = bs_dev
        self.bs_dev_plus_minus = self._get_deviation_from_ideal_value(inputarray=self.array_bvs, ideal=ideal_bs)
        dev_array, self.frequency = self._get_frequency_of_values(arraydeviations=bs_dev, stepsize=0.01)
        dev_array_plus_minus, self.frequency_plus_minus = self._get_frequency_bs_for_plus_minus_deviation(
            self.bs_dev_plus_minus, stepsize=self.stepsize_histogram)
        self.arraydev_share = self._list_to_np_array_and_divide_by_value(np.array(dev_array), ideal_bs)
        self.relativefrequency = self._list_to_np_array_and_divide_by_value(np.array(self.frequency),
                                                                            float(len(array_bvs)))

        self.arraydev_share_plus_minus = self._list_to_np_array_and_divide_by_value(np.array(dev_array_plus_minus),
                                                                                    ideal_bs)

        self.relativefrequency_plus_minus = self._list_to_np_array_and_divide_by_value(
            np.array(self.frequency_plus_minus),
            float(len(array_bvs)))
        # TODO: iclude new functions to calculate positive and neg deviation + histogram
        # TODO: safe this
        # TODO: transform it into a bar chart
        # TODO: bin this correctly +- stepsize etc


class Pauling3OverAllAnalysis(OverAllAnalysis):
    """
    Class to analyse the third Pauling rule
    """

    def run(self, show_plot=True, start_from_connections=False, save_connections=True,
            connections_folder='AnalysisConnections', start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analysis=False, save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json', start_material=None, stop_material=None, maxCN=None,
            EdgesAsAdditionalExceptions=False):
        """

        :param show_plot: if True, the pie plot, and the dependencies of the rule fulfillment on atomic radii and mean CN are shown
        :param start_from_connections: starts from saved connection files (recommended if they exist)
        :param save_connections: if true, will save the connections if not already present
        :param connections_folder: the folder, in which the connection files are saved
        :param start_from_results: if true, will only read and plot the results
        :param save_result_data: if true, will save complete result data
        :param restart_from_saved_structure_analysis: if true, will restart the structure analysis from before
        :param save_structure_analysis: if true, will save the structure analysis
        :param path_to_save: where shall the results be saved (structure analysis will be saved in an adapted path in the same directory)
        :param start_material: number of the material at which the analysis shall be started
        :param stop_material: number of the material before which the analysis shall be stopped
        :param maxCN: only polyhedra with this coordination number or lower will be considered in the analysis
        :param EdgesAsAdditionalExceptions: Edges are also treated as exceptions in the elementwise analysis
        :return:
        """

        self.start_material = start_material
        self.stop_material = stop_material
        self.connections_folder = connections_folder
        self.start_from_connections = start_from_connections
        self.save_connections = save_connections
        self.maxCN = maxCN
        self.EdgesAsAdditionalExceptions = EdgesAsAdditionalExceptions
        if not start_from_results:
            self._new_setup()
            # warum funktioniert das nicht?
            self.Plot_PSE_numbers = get_mean_CN_from_frequencies(self.All_Details)
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])
            self.connections = inputdict['connections']
            self.Plot_PSE_numbers = inputdict['PSE_Dict2']
            self.additional_info = inputdict['additional_info']

        if show_plot:
            # do the other plot here
            # print(self.Plot_PSE_numbers)
            plot1 = self._plot_influence_mean_CN(self.Plot_PSE_DICT, self.Plot_PSE_numbers,
                                                 self.EdgesAsAdditionalExceptions)
            plot1.show()
            plot2 = self._plot_influence_atomic_radii(self.Plot_PSE_DICT, self.EdgesAsAdditionalExceptions)
            plot2.show()

            plot = self._pieplot_connections(self.connections['corner'], self.connections['edge'],
                                             self.connections['face'], 'Connected Pairs of Polyhedra')
            plot.show()

            # TODO: plot that shows number of connections via corners and edges vs. atomic radii, and number of connections via faces vs. atomic radii: get radii first
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
            outputdict['PSE_Dict2'] = self.Plot_PSE_numbers
            outputdict['additional_info'] = self.additional_info
            self._save_results_to_file(outputdict, path_to_save)

        print(str(len(self.structures_fulfillingrule) / (len(self.structures_fulfillingrule) + len(
            self.structures_exceptions))) + " of all tested structures fulfill the rule.")

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling

            if save_structure_analysis:
                # self._print_to_screen_similar_structures(dict_similarstructures_fulfilling)
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections of Cations with Valences')

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections of Cations with Valences')

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections of Cations with Valences')

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections of Cations with Valences')

    def _plot_influence_atomic_radii(self, Plot_PSE_Dict: dict, EdgeAsAdditionalException=False) -> plt:
        """
        will plot number of connected polyhedra vs. atomic radii
        :param Plot_PSE_Dict: a dict with information on the rule fulfillment for each cation
        :param EdgeAsAdditonalExceptions: Edges will also be treated as additional exceptions
        :return:
        """
        # oder histogramm
        corneredge_histo = []
        face_histo = []
        for key, item in Plot_PSE_Dict.items():
            for number in range(0, item[0]):
                corneredge_histo.append(Element(key).atomic_radius)
            for number in range(0, item[1]):
                face_histo.append(Element(key).atomic_radius)

        # TODO: make a statistical test and test for the different distributions
        # TODO: utilize a smilar plot for 4th rule
        if not EdgeAsAdditionalException:
            print("Mean CN for Corners and Edge Connected Polyhedra")
            print(np.mean(corneredge_histo))
            print("Standard Deviation CN for Corners and Edge Connected Polyhedra")
            print(np.std(corneredge_histo, ddof=1))
            print("Mean CN for Face Connected Polyhedra")
            print(np.mean(face_histo))
            print("Standard Deviation CN for Face Connected Polyhedra")
            print(np.std(face_histo, ddof=1))
        else:
            print("Mean CN for Corners Connected Polyhedra")
            print(np.mean(corneredge_histo))
            print("Standard Deviation CN for Corners Connected Polyhedra")
            print(np.std(corneredge_histo, ddof=1))
            print("Mean CN for Edge and Face Connected Polyhedra")
            print(np.mean(face_histo))
            print("Standard Deviation CN for Edge and Face Connected Polyhedra")
            print(np.std(face_histo, ddof=1))

        # range noch korrekt anpassen
        plt.subplot(2, 1, 1)
        plt.hist(corneredge_histo, bins=[i / float(1000) for i in
                                         range(int(min([min(corneredge_histo), min(face_histo)]) * 1000),
                                               int(max([max(corneredge_histo), max(face_histo)]) * 1000) + 20, 20)],
                 align='mid', alpha=1)
        plt.axvline(np.mean(corneredge_histo), color='r')
        if not EdgeAsAdditionalException:
            plt.ylabel("Connections via corners and edges")
        else:
            plt.ylabel("Connections via corners")
        # plt.plot([0, 1000], [np.mean(corneredge_histo), np.mean(corneredge_histo)], 'b-')
        # n,bins,patches=plt.hist(x=corneredge_histo,bins=len(range(int(min(corneredge_histo)*100),int(max(corneredge_histo)*100)+5)))
        # y =sp.stats.norm.pdf(bins,np.mean(corneredge_histo),np.std(corneredge_histo,ddof=1))
        # plt.plot(bins,y,'r--',linewidth=1)
        #

        plt.subplot(2, 1, 2)
        plt.hist(face_histo, bins=[i / float(1000) for i in
                                   range(int(min([min(corneredge_histo), min(face_histo)]) * 1000),
                                         int(max([max(corneredge_histo), max(face_histo)]) * 1000) + 20, 20)],
                 align='mid',
                 alpha=1)
        plt.axvline(np.mean(face_histo), color='r')
        if not EdgeAsAdditionalException:
            plt.ylabel("Connections via faces")
        else:
            plt.ylabel("Connections via edges and faces")

        plt.xlabel("Atomic radius in Angstrom")
        # plt.plot([0,1000],[np.mean(face_histo),np.mean(face_histo)],'b-')
        return plt

    def _plot_influence_mean_CN(self, Plot_PSE_Dict: dict, Plot_PSE_CN_Info: dict,
                                EdgeAsAdditionalException=False) -> plt:
        """
        will plot number of conncted polyhedra vs. mean CN
        :param Plot_PSE_Dict: dict with info on the rule fulfillment for each cation
        :param Plot_PSE_CN_Info: dict with info for each cation on the mean CN
        :param EdgeAsAdditonalExceptions: Edges will also be treated as additional exceptions
        :return: plot
        """
        # oder histogramm
        corneredge_histo = []
        face_histo = []
        for key, item in Plot_PSE_Dict.items():
            for number in range(0, item[0]):
                corneredge_histo.append(Plot_PSE_CN_Info[key])
            for number in range(0, item[1]):
                face_histo.append(Plot_PSE_CN_Info[key])

        # TODO: make a statistical test and test for the different distributions
        # TODO: utilize a smilar plot for 4th rule
        print(np.mean(corneredge_histo))
        print(np.std(corneredge_histo, ddof=1))
        print(np.mean(face_histo))
        print(np.std(face_histo, ddof=1))

        # range noch korrekt anpassen
        plt.subplot(2, 1, 1)
        plt.hist(corneredge_histo, bins=[i / float(1000) for i in
                                         range(int(min([min(corneredge_histo), min(face_histo)]) * 1000),
                                               int(max([max(corneredge_histo), max(face_histo)]) * 1000) + 20, 20)],
                 align='mid', alpha=1)
        plt.axvline(np.mean(corneredge_histo), color='r')
        if not EdgeAsAdditionalException:
            plt.ylabel("Connections via corners and edges")
        else:
            plt.ylabel("Connections via corners")
        # plt.plot([0, 1000], [np.mean(corneredge_histo), np.mean(corneredge_histo)], 'b-')
        # n,bins,patches=plt.hist(x=corneredge_histo,bins=len(range(int(min(corneredge_histo)*100),int(max(corneredge_histo)*100)+5)))
        # y =sp.stats.norm.pdf(bins,np.mean(corneredge_histo),np.std(corneredge_histo,ddof=1))
        # plt.plot(bins,y,'r--',linewidth=1)
        #

        plt.subplot(2, 1, 2)
        plt.hist(face_histo, bins=[i / float(1000) for i in
                                   range(int(min([min(corneredge_histo), min(face_histo)]) * 1000),
                                         int(max([max(corneredge_histo), max(face_histo)]) * 1000) + 20, 20)],
                 align='mid',
                 alpha=1)
        plt.axvline(np.mean(face_histo), color='r')
        if not EdgeAsAdditionalException:
            plt.ylabel("Connections via faces")
        else:
            plt.ylabel("Connections via edges and faces")
        plt.xlabel("Mean of the CN")
        # plt.plot([0,1000],[np.mean(face_histo),np.mean(face_histo)],'b-')
        return plt

    def _new_setup(self):
        """
        will setup calculations from scratch
        :return:
        """
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.additional_info = {}
        self.Plot_PSE_DICT = {}
        self.connections = {"corner": 0, "edge": 0, "face": 0}
        self.present_env = {}
        self.All_Details = {}

        for mat in list_mat:
            # print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)
            pauling0 = Pauling0(lse)
            pauling1_limit = FrequencyEnvironmentPauling1(lse=lse)
            pauling3 = Pauling3()

            if not self.start_from_connections or not os.path.isfile(
                    os.path.join(self.connections_folder, mat + '.json')):
                pauling3.newsetup(lse, filename=mat + '.json', save_to_file=self.save_connections,
                                  foldername=self.connections_folder, distance=8.0)
            else:
                pauling3.from_file(filename=mat + '.json', foldername=self.connections_folder)
            try:
                if pauling3.is_fulfilled(maximumCN=self.maxCN):
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)

            Details = pauling3.get_details(maximumCN=self.maxCN)
            self.additional_info[mat] = Details["species"]

            New_Details = self._reformat_details(Details['species'],
                                                 EdgesAsAdditionalExceptions=self.EdgesAsAdditionalExceptions)
            self._add_dict_cat_dependency(self.Plot_PSE_DICT, New_Details)

            self._sum_connections(self.connections, Details)

            self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                          number_of_elements_to_add=1)

            Details2 = pauling1_limit.get_details()

            self._add_dict_cat_dependency(start_dict=self.All_Details, dict_to_add=Details2,
                                          number_of_elements_to_add=4)

    def _reformat_details(self, Details: dict, EdgesAsAdditionalExceptions=False) -> dict:
        """
        #TODO: test this again!
        reformats dicts
        :param Details:
        :return:
        """

        New_Details = {}
        for key, item in Details.items():
            if not key in New_Details:
                New_Details[key] = [0, 0]
            for item2 in item.values():
                New_Details[key][0] += item2['corner']
                if not EdgesAsAdditionalExceptions:
                    New_Details[key][0] += item2['edge']
                else:
                    New_Details[key][1] += item2['edge']
                New_Details[key][1] += item2['face']

        return New_Details

    def _sum_connections(self, start_dict: dict, Details: dict):
        """
        sums dicts
        :param start_dict: here, the other values will be added
        :param Details: the dict that will be added
        :return:
        """
        start_dict["corner"] += Details["corner"]
        start_dict["edge"] += Details["edge"]
        start_dict["face"] += Details["face"]


class Pauling4OverAllAnalysis(OverAllAnalysis):
    """
    Class to analyse the fourth rule
    """

    def run(self, show_plot=True, start_from_connections=False, save_connections=True,
            connections_folder='AnalysisConnections', start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analysis=False, save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json', start_material=None, stop_material=None):
        """

        :param show_plot: if True, main plots for fourth rule will be shown
        :param start_from_connections: if true, will start from the calculated connections
        :param save_connections: if true, will save the connections
        :param connections_folder: folder, where the connection files are saved
        :param start_from_results: if true, will start from result file
        :param save_result_data: if true, will save result data
        :param restart_from_saved_structure_analysis: if true, will restart from a previous structure analysis
        :param save_structure_analysis: if true, will save the structure analysis in the same folder as the other results
        :param path_to_save: path where the files shall be saved
        :param start_material: number of the material at which the analysis will be started
        :param stop_material: number of the material before which the analysis will be stopped
        :return:
        """

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
            self.additional_info = inputdict['additional_info']

        print(
            str(len(self.structures_fulfillingrule) + len(self.structures_exceptions)) + " structures were considered.")
        # print(str(len(self.structures_cannot_be_evaluated))+ ' were not considered')
        if show_plot:
            plot = self._fourthrule_plot(self.Dict_val, option='val', vmin=0.7, vmax=1.0)
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
            outputdict['additional_info'] = self.additional_info
            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling

            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections depending on element, CN, and valence')

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections depending on element, CN, and valence')

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections depending on element, CN, and valence')

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv",
                                                       add_info=self.additional_info,
                                                       name_add_info='Connections depending on element, CN, and valence')

    def _new_setup(self):
        """
        will start analysis from scratch
        :return:
        """
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries,
                                            start_material=self.start_material, stop_material=self.stop_material)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.Dict_val = {}
        self.Dict_CN = {}
        self.Plot_PSE_DICT = {}
        self.present_env = {}
        self.additional_info = {}
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
                self.additional_info[mat] = Details['elementwise']

                # think about symmetry - there might be something wrong here
                New_Details_val = self._reformat_details_val(Details)
                New_Details_CN = self._reformat_details_CN(Details)
                # print(New_Details_val)
                # print(New_Details_CN)
                # # tests if matrices are symmetric!
                # if not self._check_details_symmetric(New_Details_val):
                #     raise ValueError
                # if not self._check_details_symmetric(New_Details_CN):
                #     raise ValueError

                self._add_dict_CN_val(self.Dict_val, New_Details_val)
                self._add_dict_CN_val(self.Dict_CN, New_Details_CN)

                Elementwise_analysis = self._reformat_details_elementwise(Details)
                # print(Elementwise_analysis)
                self._add_dict_cat_dependency(self.Plot_PSE_DICT, Elementwise_analysis)
                # print(Details['elementwise'])
                # print(Details['maxval'])
                # print(Details['minCN'])
                # print(Elementwise_analysis)

                self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                              number_of_elements_to_add=1)

            except RuleCannotBeAnalyzedError:
                pass

    def _reformat_details_elementwise(self, Details: dict) -> dict:
        """
        will reformat a dict
        :param Details: dict
        :return: dict
        """
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

    def _add_dict_CN_val(self, start_dict: dict, to_add_dict: dict):
        """
        will sum two dicts
        :param start_dict: dict to which the new dict will be added
        :param to_add_dict: dict to add
        :return:
        """
        for key, item in to_add_dict.items():
            if key not in start_dict:
                start_dict[key] = {}
            for key2, item2 in item.items():
                if key2 not in start_dict[key]:
                    start_dict[key][key2] = [0, 0]
                start_dict[key][key2][0] += to_add_dict[key][key2][0]
                start_dict[key][key2][1] += to_add_dict[key][key2][1]

    def _reformat_details_val(self, Details: dict) -> dict:
        """
        reformats dict
        :param Details:  dict
        :return: dict
        """
        # do a simple copy
        New = {}
        for key1, item1 in Details.items():
            if "val1" in key1:
                for key2, item2 in item1.items():
                    for key3, item3 in item2.items():
                        for key4, item4 in item3.items():
                            if not key4 in New:
                                New[key4] = {}
                            if not key3 in New[key4]:
                                New[key4][key3] = {}
                            if not key2 in New[key4][key3]:
                                New[key4][key3][key2] = {}
                            if not key1 in New[key4][key3][key2]:
                                New[key4][key3][key2][key1] = item4

        # similar to CN method
        New_Details = {}
        for key1, item1 in New.items():
            CN1 = key1.split(":")[1]
            for key2, item2 in item1.items():

                CN2 = key2.split(":")[1]
                for key3, item3 in item2.items():
                    val1 = key3.split(":")[1]
                    if val1 not in New_Details:
                        New_Details[val1] = {}

                    for key4, item4 in item3.items():
                        val2 = key4.split(":")[1]
                        if val2 not in New_Details[val1]:
                            New_Details[val1][val2] = [0, 0]
                        if val2 not in New_Details:
                            New_Details[val2] = {}
                        if val1 not in New_Details[val2]:
                            New_Details[val2][val1] = [0, 0]

                        New_Details[val1][val2][0] += item4["no"]
                        New_Details[val1][val2][1] += item4["corner"]
                        New_Details[val1][val2][1] += item4["edge"]
                        New_Details[val1][val2][1] += item4["face"]
                        # to arrive at a correct sum of the connections
                        # typically this will result in the number of atoms that are connected, not in the number of connections
                        if not (val1 == val2 and CN1 != CN2):
                            New_Details[val2][val1][0] += item4["no"]
                            New_Details[val2][val1][1] += item4["corner"]
                            New_Details[val2][val1][1] += item4["edge"]
                            New_Details[val2][val1][1] += item4["face"]
        return New_Details

    def _reformat_details_CN(self, Details: dict) -> dict:
        """
        will reformat the dict
        :param Details: dict
        :return:  dict
        """
        New_Details = {}

        for key1, item1 in Details.items():

            if "val1" in key1:
                val1 = key1.split(":")[1]
                for key2, item2 in item1.items():
                    val2 = key2.split(":")[1]
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
                            # very similar to val method
                            # make sure the number of connections is correct
                            # also the dict has to be symmetric under exchange of CN
                            if not (CN1 == CN2 and val1 != val2):
                                New_Details[CN2][CN1][0] += item4["no"]
                                New_Details[CN2][CN1][1] += item4["corner"]
                                New_Details[CN2][CN1][1] += item4["edge"]
                                New_Details[CN2][CN1][1] += item4["face"]

        return New_Details

    def _plot_get_max(self, inputdict: dict) -> float:
        """
        get max  value of keys of dict
        :param inputdict
        :return: float
        """
        return max([int(x) for x in inputdict.keys()])

    def _plot_init_np(self, maxvalue: int) -> np.full:
        """
        initalizes numpy.full
        :param maxvalue: int that determines how large np.full will be
        :return: numpy.full
        """
        return np.full((maxvalue, maxvalue), np.nan)

    def _plot_fill_np(self, PlotEndDict, inputdict, lowervalue):
        """
        will fill in the entries of the plot dict that will then be plotted
        :param PlotEndDict: dict that will later be plotted
        :param inputdict: input dictionary
        :param lowervalue: only if there are more than lowervalue entries, there will be an entry in the ouptut dict
        :return:
        """
        for key, item in inputdict.items():
            for key2, item2 in item.items():
                if item2[0] + item2[1] > lowervalue:
                    PlotEndDict[int(key)][int(key2)] = float(item2[0]) / float(item2[0] + item2[1])

    def _fourthrule_plot(self, inputdict, option='CN', vmin=0.0, vmax=1.0, lowervalue=50, leaveoutCN13=True):
        """
        #leaves out CN 13 for option "CN" -> too few entries
        :param inputdict: dict that will be plotted
        :param option: 'CN' or 'val' -> changes labels of axis dynamically
        :param vmin: lower limit of plotted z values
        :param vmax: upper limit of plotted z values
        :param: lowervalue: only if more than lowervalue polyhedra were considered, it will be plotted
        :param: leaveoutCN13: if  true, CN=13 will be left out from plot
        :return:
        """
        # get valmax
        valCNmax = self._plot_get_max(inputdict)

        if option == 'val':
            PlotEndDict = self._plot_init_np(valCNmax + 1)
        elif option == 'CN':
            if not leaveoutCN13:
                PlotEndDict = self._plot_init_np(valCNmax + 2)
            else:
                PlotEndDict = self._plot_init_np(valCNmax + 2)

        self._plot_fill_np(PlotEndDict, inputdict, lowervalue)

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
    """
    Class to analyse the 5th rule
    """

    def run(self, remove_elements_low_entropy=True, show_plot=True, start_from_connections=False, save_connections=True,
            connections_folder='AnalysisConnections_5thRule', start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analysis=False, save_structure_analysis=True,
            path_to_save='Results/Results_Second_Rule.json',
            threshold_remove_elements=0.95, start_material=None, stop_material=None):
        """
        will run the overall analysis
        :param remove_elements_low_entropy: removes elements that have a low shannon entropy
        :param show_plot: if true, will display all relevant plots
        :param start_from_connections: if true, will start from connection files
        :param save_connections: if true, will save connections in files
        :param connections_folder: folder, in which the connections are saved
        :param start_from_results: if true, will start from saved results
        :param save_result_data: if true, will save results
        :param restart_from_saved_structure_analysis: will restart from a structure analysis that was begun
        :param save_structure_analysis: will save structure analysis, if true
        :param path_to_save: where will the results be saved
        :param threshold_remove_elements: will determine which elements will be removed according to shannon entropy
        :param start_material: number, from which the analysis will be started
        :param stop_material: number, before whch the analysis will be stopped
        :return:
        """
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
                                           plot_element_dependend_analysis=False,
                                           list_of_materials_to_investigate=self.list_of_materials_to_investigate)
                newclass.run(start_from_results=False, save_result_data=False)
                self.list_to_remove = [el for el, item in newclass.Plot_PSE_entropy.items() if
                                       item > threshold_remove_elements]
            else:
                self.list_to_remove = []
            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.Plot_PSE_DICT = inputdict['PSE_Dict']
            self.present_env = Counter(inputdict['Counter_cation'])
            self.list_to_remove = inputdict['elements_low_entropy']

        print(
            str(len(self.structures_fulfillingrule) + len(self.structures_exceptions)) + " structures were considered.")
        # print(str(len(self.structures_cannot_be_evaluated)) + ' were not considered')

        if show_plot:
            # print(len(self.structures_fulfillingrule))
            # print(len(self.structures_exceptions))
            plot = self._fifth_rule_plot(okay=len(self.structures_fulfillingrule),
                                         notokay=len(self.structures_exceptions))
            plot.show()

        if not start_from_results:
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
            outputdict['elements_low_entropy'] = self.list_to_remove
            self._save_results_to_file(outputdict, path_to_save)

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            self.dict_similarstructures_exceptions = dict_similarstructures_exceptions
            self.dict_similarstructures_fulfilling = dict_similarstructures_fulfilling

            if save_structure_analysis:
                self._print_to_file_similar_structures(dict_similarstructures_exceptions,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling,
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.yaml")

                self._print_to_file_similar_structures(dict_similarstructures_exceptions, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structural_exceptions_readable.csv")

                self._print_to_file_similar_structures(dict_similarstructures_fulfilling, fmt='csv',
                                                       filename=path_to_save.split('.')[
                                                                    0] + "_structures_fulfilling_readable.csv")

    def _new_setup(self):
        """
        will start the analysis from scratch
        :return:
        """
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
            # print(mat)
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
                # print(Details)
                New_Details = self._reformat_details_elementdependency(Details)

                self._add_dict_cat_dependency(self.Plot_PSE_DICT, New_Details)
                self._add_dict_cat_dependency(self.present_env, pauling0.get_cations_in_structure(),
                                              number_of_elements_to_add=1)


            except RuleCannotBeAnalyzedError:
                pass

            # else:
            #     counter_not_primitive+=1
            #     print(counter_not_primitive)

    def _fifth_rule_plot(self, okay: int, notokay: int, title="Tested Structures") -> plt:
        """
        will plot a pie plot
        :param okay: number of structures fulfilling the rule
        :param notokay: number of structures that do not fulfill the rule
        :param title: title of the plot
        :return: will return a plot
        """
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

    def _reformat_details_elementdependency(self, details: dict) -> dict:
        """
        will reformat a dict to make the analysis easier
        :param details: dict
        :return: dict
        """
        new_dict = {}
        for key in details:
            if not key in new_dict:
                new_dict[key] = [details[key]['fulfilled'], details[key]['not_fulfilled']]
        return new_dict


class AllPaulingOverAllAnalysis(OverAllAnalysis):
    """
    Class to analyse all 5 Pauling rules
    """

    # TODO: overwrite init to avoid confusions

    def run(self, remove_elements_low_entropy=False, start_from_connections=False,
            save_connections=True, connections_folder34='AnalysisConnections',
            connections_folder5='AnalysisConnections_5thRule',
            start_from_results=False, save_result_data=True,
            restart_from_saved_structure_analysis=False, save_structure_analysis=True,
            path_to_save='Results/Results_AllRules.json', threshold_remove_elements=0.95, start_material=None,
            stop_material=None, adapt_first_fourth_and_fifth_rules=False, ignore_first_rule=True,
            ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False, ignore_fifth_rule=False,
            remove_structures_with_CN_larger_8=False, compute_elementwise=False):
        """
        will run the overall analysis
        :param remove_elements_low_entropy: will remove the elements with low shannon entropy
        :param start_from_connections: if true, will start from connections files
        :param save_connections: if true, will save the connections files
        :param connections_folder34: folder for 3rd and 4th rule files
        :param connections_folder5: folder for 5th rule files
        :param start_from_results: if true, will start from results
        :param save_result_data: if true, will save results
        :param restart_from_saved_structure_analysis: if true, will restart from a started structure analysis
        :param save_structure_analysis: if true, will save the structure analysis
        :param path_to_save: path, where results will be saved
        :param threshold_remove_elements: will decide, which elements will be removed according to shannon entropy
        :param start_material: material number, at which the analysis will be started
        :param stop_material: material number, before which the analysis will be stated
        :param adapt_first_fourth_and_fifth_rules: will adapt the first, fourth, fifth rule (different from before)
        :param ignore_first_rule: if true, will skip analysis of first rule
        :param ignore_second_rule: if true, will skip analysis of second rule
        :param ignore_third_rule: if true, will skip analysis of third rule
        :param ignore_fourth_rule: if true, will skip analysis of fourth rule
        :param ignore_fifth_rule: if true, will skip analysis of fifth rule
        :param remove_structures_with_CN_larger_8: if true, will remove all structures with CN>8
        :param will compute elementwise fulfillment of the rule even if plotting is switched off
        :return:
        """

        self.start_material = start_material
        self.stop_material = stop_material
        self.ignore_first_rule_and_test_criteria_for_rule_four_and_five = adapt_first_fourth_and_fifth_rules
        self.ignore_first_rule = ignore_first_rule
        self.ignore_second_rule = ignore_second_rule
        self.ignore_third_rule = ignore_third_rule
        self.ignore_fourth_rule = ignore_fourth_rule
        self.ignore_fifth_rule = ignore_fifth_rule

        self.remove_elements_low_entropy = remove_elements_low_entropy
        self.connections_folder34 = connections_folder34
        self.connections_folder5 = connections_folder5
        self.start_from_connections = start_from_connections
        self.save_connections = save_connections
        self.threshold_remove_elements = threshold_remove_elements

        if not start_from_results:
            if self.remove_elements_low_entropy:
                newclass = Pauling1Entropy(source=self.source, onlybinaries=self.onlybinaries,
                                           plot_element_dependend_analysis=False,
                                           list_of_materials_to_investigate=self.list_of_materials_to_investigate,
                                           start_material=self.start_material, stop_material=self.stop_material)
                newclass.run(start_from_results=False, save_result_data=False)
                self.list_to_remove = [el for el, item in newclass.Plot_PSE_entropy.items() if
                                       item > self.threshold_remove_elements]
            else:
                self.list_to_remove = []

            self._new_setup()
        else:
            inputdict = self._get_precomputed_results(path_to_save)
            self.structures_fulfillingrule = inputdict['fullingstructures']
            self.structures_exceptions = inputdict['exceptions']
            self.structures_cannot_be_evaluated = inputdict['nottested']
            self.list_to_remove = inputdict['list_to_remove']

        if save_result_data and not start_from_results:
            outputdict = {}
            outputdict['fullingstructures'] = self.structures_fulfillingrule
            outputdict['exceptions'] = self.structures_exceptions
            outputdict['nottested'] = self.structures_cannot_be_evaluated
            outputdict['list_to_remove'] = self.list_to_remove
            self._save_results_to_file(outputdict, path_to_save)
        # TODO: test both parts -> remove rule 3,4 , also test this for structures that are to too big
        if remove_structures_with_CN_larger_8:
            structures_fulfilling_2 = []
            for mat in self.structures_fulfillingrule:

                lse = self._get_lse_from_folder(mat=mat, source=self.source)
                CNtoobig = False
                for isite, site_envs in enumerate(lse.coordination_environments):
                    if site_envs != None:
                        if len(site_envs) > 0:
                            CN = int(site_envs[0]['ce_symbol'].split(':')[1])
                            if CN > 8:
                                CNtoobig = True
                                break
                if not CNtoobig:
                    structures_fulfilling_2.append(mat)
            self.structures_fulfillingrule = structures_fulfilling_2

            structures_exceptions_2 = []
            for mat in self.structures_exceptions:

                lse = self._get_lse_from_folder(mat=mat, source=self.source)
                CNtoobig = False
                for isite, site_envs in enumerate(lse.coordination_environments):
                    if site_envs != None:
                        if len(site_envs) > 0:
                            CN = int(site_envs[0]['ce_symbol'].split(':')[1])
                            if CN > 8:
                                CNtoobig = True
                                break
                if not CNtoobig:
                    structures_exceptions_2.append(mat)
            self.structures_exceptions = structures_exceptions_2

        if len(self.structures_exceptions) + len(self.structures_fulfillingrule) > 0:
            print("Percentage of fulfillment")
            print(len(self.structures_fulfillingrule) / (
                    len(self.structures_exceptions) + len(self.structures_fulfillingrule)))

            self.percentage_structures_fulfilling = len(self.structures_fulfillingrule) / (
                    len(self.structures_exceptions) + len(self.structures_fulfillingrule))
            print("Exceptions")
            print(len(self.structures_exceptions))
            print("Fulfilling")
            print(len(self.structures_fulfillingrule))

        print(
            str(len(self.structures_fulfillingrule) + len(self.structures_exceptions)) + " structures were considered.")
        # print(str(len(self.structures_cannot_be_evaluated)) + ' were not considered')

        self.Plot_PSE = {}
        # Analyze the compositions and compute the fulfillment per element here
        # TODO: test this part in test code
        if self.plot_element_dependend_analysis or compute_elementwise:
            for mat in self.structures_fulfillingrule:
                lse = self._get_lse_from_folder(mat, source=self.source)
                elements = lse.structure.composition.elements
                # print(elements)
                for element in elements:
                    if str(element) != 'O':
                        if not str(element) in self.Plot_PSE:
                            self.Plot_PSE[str(element)] = [0, 0]
                        self.Plot_PSE[str(element)][0] += 1
            for mat in self.structures_exceptions:
                lse = self._get_lse_from_folder(mat, source=self.source)
                elements = lse.structure.composition.elements
                # print(elements)
                for element in elements:
                    if str(element) != 'O':
                        if not str(element) in self.Plot_PSE:
                            self.Plot_PSE[str(element)] = [0, 0]
                        self.Plot_PSE[str(element)][1] += 1
        if self.plot_element_dependend_analysis:
            plot = self._plot_PSE(self.Plot_PSE,
                                  lowest_number_of_environments_considered=self.lowest_number_environments_for_plot,
                                  lowerlimit=self.lower_limit_plot, upperlimit=self.upper_limit_plot)
            plot.show()

        if self.analyse_structures:
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source,
                                                                             save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structures_fulfilling.json",
                                                                             fetch_results_only=restart_from_saved_structure_analysis,
                                                                             start_from_Matching=self.use_prematching)

        if save_structure_analysis and self.analyse_structures:
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
        """
        will start overall analysis from scratch
        :return:
        """
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
        # print(self.list_to_remove)
        for mat in list_mat:
            lse = self._get_lse_from_folder(mat, source=self.source)
            pauling1 = Pauling1(lse, filenameradii='../Assessment/Should_not_be_changed/univalent_cat_radii.json',
                                onlylowerlimit=False)
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
            # TODO: leave rule  1 in test
            if not self.ignore_first_rule_and_test_criteria_for_rule_four_and_five:
                try:
                    if not self.ignore_first_rule:
                        result1 = pauling1.is_fulfilled()
                    else:
                        result1 = True
                    if not self.ignore_second_rule:
                        result2 = pauling2.is_fulfilled()
                    else:
                        result2 = True
                    if not self.ignore_third_rule:
                        result3 = pauling3.is_fulfilled()
                    else:
                        result3 = True
                    if not self.ignore_fourth_rule:
                        result4 = pauling4.is_fulfilled()
                    else:
                        result4 = True
                    if not self.ignore_fifth_rule:
                        result5 = pauling5.is_fulfilled(leave_out_list=self.list_to_remove)
                    else:
                        result5 = True

                    if result1 and result2 and result3 and result4 and result5:
                        self.structures_fulfillingrule.append(mat)
                    else:
                        self.structures_exceptions.append(mat)
                except RuleCannotBeAnalyzedError:
                    self.structures_cannot_be_evaluated.append(mat)
            else:
                try:

                    # TODO: include more options

                    if not self.ignore_first_rule:
                        Details_Pauling1 = pauling1.get_details()
                        if Details_Pauling1['Env_notfulfilled'] == 0 and Details_Pauling1['Env_fulfilled'] > 0:
                            result1 = True
                        elif Details_Pauling1['Env_notfulfilled'] == 0 and Details_Pauling1["Env_fulfilled"] == 0:
                            raise RuleCannotBeAnalyzedError
                        else:
                            result1 = False
                    else:
                        result1 = True

                    # TODO: include removing of rules!
                    if not self.ignore_second_rule:
                        result2 = pauling2.is_fulfilled()
                    else:
                        result2 = True

                    if not self.ignore_third_rule:
                        result3 = pauling3.is_fulfilled()
                    else:
                        result3 = True

                    if not self.ignore_fourth_rule:
                        try:
                            result4 = pauling4.is_fulfilled()
                        except RuleCannotBeAnalyzedError:
                            result4 = True
                    else:
                        result4 = True

                    if not self.ignore_fifth_rule:
                        try:
                            result5 = pauling5.is_fulfilled()
                        except RuleCannotBeAnalyzedError:
                            result5 = True
                    else:
                        result5 = True

                    if result1 and result2 and result3 and result4 and result5:
                        self.structures_fulfillingrule.append(mat)
                    else:
                        self.structures_exceptions.append(mat)

                except RuleCannotBeAnalyzedError:
                    self.structures_cannot_be_evaluated.append(mat)


class AllPaulingOverAllAnalysis_Final_Summary(OverAllAnalysis):
    """
    Class to analyse 4 rules and the dependency of each rule on the rule fulfillment
    """

    # TODO: overwrite init to avoid confusions

    def run(self, remove_elements_low_entropy=False, start_from_connections=False,
            save_connections=True, connections_folder34='AnalysisConnections',
            connections_folder5='AnalysisConnections_5thRule',
            start_from_results=False, save_result_data=True,
            path_to_save='Results/Results_AllRules_Final_plot.json', threshold_remove_elements=0.95,
            start_material=None,
            stop_material=None, plot_result=True):
        """

        :param remove_elements_low_entropy: elements with low shannon entropy will be removed
        :param start_from_connections: if true, will start from connections files
        :param save_connections: if true, save connection files
        :param connections_folder34: folder for 3rd and 4th rule files
        :param connections_folder5: folder for 5th rule files
        :param start_from_results: if true, start from results
        :param save_result_data: if true, saves results
        :param path_to_save: path to save the files
        :param threshold_remove_elements: decides which elements will be removed
        :param start_material: material number, at which the analysis is started
        :param stop_material: material number, beore the analysis will be stopped
        :param plot_result: if true, will plot results
        :return:
        """
        # general parameters, could be used in run instead
        adapt_first_fourth_and_fifth_rules = True
        ignore_first_rule = True

        # TODO: lad den kram iwie vom file
        if start_from_results:
            inputdict = self._get_precomputed_results(path_to_save)
            self.means_CN_all = inputdict["means_CN_all"]
            self.means_CN_smaller9 = inputdict["means_CN_smaller9"]
        else:

            self.means_CN_all = []
            self.means_CN_smaller9 = []
            newclass = AllPaulingOverAllAnalysis(source=self.source, onlybinaries=self.onlybinaries,
                                                 analyse_structures=False, use_prematching=True,
                                                 plot_element_dependend_analysis=self.plot_element_dependend_analysis,
                                                 list_of_materials_to_investigate=self.list_of_materials_to_investigate)

            # all rules
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=False)

            self.means_CN_all.append(newclass.percentage_structures_fulfilling)

            # leave out second rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=True, ignore_third_rule=False, ignore_fourth_rule=False,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=False)

            self.means_CN_all.append(newclass.percentage_structures_fulfilling)

            # leave out third rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=True, ignore_fourth_rule=False,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=False)

            self.means_CN_all.append(newclass.percentage_structures_fulfilling)

            # leave out fourth rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=True,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=False)

            self.means_CN_all.append(newclass.percentage_structures_fulfilling)

            # leave out fifth rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False,
                         ignore_fifth_rule=True, remove_structures_with_CN_larger_8=False)

            self.means_CN_all.append(newclass.percentage_structures_fulfilling)

            # CN<=8
            # all rules
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=True)

            self.means_CN_smaller9.append(newclass.percentage_structures_fulfilling)

            # leave out second rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=True, ignore_third_rule=False, ignore_fourth_rule=False,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=True)

            self.means_CN_smaller9.append(newclass.percentage_structures_fulfilling)

            # leave out third rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=True, ignore_fourth_rule=False,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=True)

            self.means_CN_smaller9.append(newclass.percentage_structures_fulfilling)

            # leave out fourth rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=True,
                         ignore_fifth_rule=False, remove_structures_with_CN_larger_8=True)

            self.means_CN_smaller9.append(newclass.percentage_structures_fulfilling)

            # leave out fifth rule
            newclass.run(remove_elements_low_entropy=remove_elements_low_entropy,
                         start_from_connections=start_from_connections, save_connections=save_connections,
                         connections_folder34=connections_folder34, connections_folder5=connections_folder5,
                         start_from_results=False, save_result_data=False,
                         restart_from_saved_structure_analysis=False, save_structure_analysis=False,
                         path_to_save='', start_material=start_material,
                         stop_material=stop_material,
                         threshold_remove_elements=threshold_remove_elements,
                         adapt_first_fourth_and_fifth_rules=adapt_first_fourth_and_fifth_rules,
                         ignore_first_rule=ignore_first_rule,
                         ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False,
                         ignore_fifth_rule=True, remove_structures_with_CN_larger_8=True)

            self.means_CN_smaller9.append(newclass.percentage_structures_fulfilling)

        if save_result_data and not start_from_results:
            outputdict = {}
            #     outputdict['fullingstructures'] = self.structures_fulfillingrule
            #     outputdict['exceptions'] = self.structures_exceptions
            #     outputdict['nottested'] = self.structures_cannot_be_evaluated
            #     outputdict['list_to_remove'] = self.list_to_remove
            outputdict["means_CN_all"] = self.means_CN_all
            outputdict["means_CN_smaller9"] = self.means_CN_smaller9
            self._save_results_to_file(outputdict, path_to_save)

        if plot_result:
            plt = self.final_plot()
            plt.show()
        # todo: speicher den kram in file

    def final_plot(self) -> plt:
        """

        :return:
        """
        import numpy as np
        import matplotlib.pyplot as plt

        # data to plot

        means_CN_all = self.means_CN_all
        means_CN_smaller8 = self.means_CN_smaller9
        n_groups = len(means_CN_all)
        # create plot
        fig, ax = plt.subplots()
        index = np.arange(n_groups)
        bar_width = 0.35
        opacity = 0.8

        rects1 = plt.bar(index, means_CN_smaller8, bar_width,
                         alpha=opacity,
                         color='b',
                         label='CN <=8')

        rects2 = plt.bar(index + bar_width, means_CN_all, bar_width,
                         alpha=opacity,
                         color='g',
                         label='all CN')

        plt.xlabel('Tested Rules')
        plt.ylabel('Percentage of Structure that fulfill rules')
        plt.title('Fulfillment of the Rules')
        plt.xticks(index + 0.5 * bar_width, ('2, 3, 4, 5', '3, 4, 5', '2, 4, 5', '2, 3, 5', '2, 3, 4'))
        plt.legend()

        plt.tight_layout()
        return plt

from PaulingRules import Pauling1, Pauling2, Pauling3, Pauling4, Pauling5, RuleCannotBeAnalyzedError
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from PlotClasses import PlotterPSE
from pymatgen.analysis.structure_matcher import StructureMatcher, FrameworkComparator
import json
import numpy as np
from collections import OrderedDict
from collections import Counter
import matplotlib
import matplotlib.pyplot as plt


class OverAllAnalysis:

    def __init__(self, source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                 lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                 analyse_structures=True):
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
        self.analyse_structures = analyse_structures
        self.lower_limit_plot = lower_limit_plot
        self.upper_limit_plot = upper_limit_plot
        # could include a function that starts plotting from saved data?
        # how should one do this in the best way

    def _get_list_materials(self, source='MP', onlybinaries=False):
        if source == 'MP':
            with open("Should_not_be_changed/allmaterials.json", "r") as f:
                list_compound_dict = json.load(f)
                list_compound = list_compound_dict["is_clear_compounds"]

        elif source == 'MP_very_symmetric':
            with open("Should_not_be_changed/ce_fraction_0.95plus_csm_0.1plus_eabovehull0.025plus_discardS1.json",
                      "r") as f:
                list_compound_dict = json.load(f)
                list_compound = list_compound_dict["is_clear_compounds"]

        if onlybinaries:
            list_compound_binaries=[]
            for mat in list_compound:
                numberofcat=0
                lse=self._get_lse_from_folder(mat=mat,source=source)
                for el in lse.structure.composition:
                    if str(el)!='O':
                        numberofcat+=1
                if numberofcat==1:
                    list_compound_binaries.append(mat)

            return list_compound_binaries

        else:
            return list_compound

    def _get_lse_from_folder(self, mat, source='MP'):
        if source == 'MP' or source == 'MP_very_symmetric':
            with open("../../DB_chem_env/1st_rule/" + mat + ".json", 'r') as f:
                data = json.load(f)

        lse = LightStructureEnvironments.from_dict(data)
        return lse

    def _plot_PSE(self, Dict_to_Plot, lowest_number_of_environments_considered, xlim=[1, 18], ylim=[1, 9], lowerlimit=0,
                  upperlimit=1, counter_cations_env=None):
        # TODO: make sure you count numbers of environments and not numbers of atoms that are counted several times

        plotterpse = PlotterPSE(valuestoplot=Dict_to_Plot, counter_cations_env=counter_cations_env)
        # atoms to consider for plot instead?

        plt = plotterpse.get_plot(xlim=xlim, ylim=ylim,
                                  lowest_number_of_environments_considered=lowest_number_of_environments_considered)
        plotterpse.set_colormap(lowerlimit, upperlimit)
        return plt

    def _save_results_to_file(self, dict_to_save, path):
        with open(path, 'w') as f:
            json.dump(dict_to_save, f)

    def _get_precomputed_results(self, path):
        with open(path, 'r') as f:
            dict_here = json.load(f)
        return dict_here

    def _add_dict_cat_dependency(self, start_dict, dict_to_add, number_of_elements_to_add=2):
        for key, item in dict_to_add.items():
            if not key in start_dict:
                if number_of_elements_to_add == 1:
                    start_dict[key] = dict_to_add[key]
                elif number_of_elements_to_add == 2:
                    start_dict[key] = dict_to_add[key].copy()
            else:
                if number_of_elements_to_add == 1:
                    start_dict[key] += dict_to_add[key]
                elif number_of_elements_to_add == 2:
                    start_dict[key][0] += dict_to_add[key][0]
                    start_dict[key][1] += dict_to_add[key][1]

    def _get_similar_structures(self, list_mat_id, source='MP', save_to_file=True,
                                path_to_save='Similar_Structures.json',fetch_results_only=False,
                                start_from_Matching=False):

        # TODO: use matching of all structures as an Input as well to speed this up
        # TODO: use further information from materials project (e.g., ICSD numbers)
        # TODO: start from file and compare with list that comes in
        # ,path_to_further_information='information_on_mat.json'
        if not start_from_Matching and not fetch_results_only:
            dictstructures = {}
            information_mat = {}
            for mat in list_mat_id:
                lse = self._get_lse_from_folder(mat, source=source)
                structure = lse.structure

                if mat not in information_mat:
                    # TODO: think which information might be interesting as well
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
                    if found:
                        dictstructures[str(foundmat)].append(mat)

                    dictstructures = OrderedDict(sorted(dictstructures.items(), key=lambda t: len(t[1]), reverse=True))

                else:
                    dictstructures[str(mat)] = []
                    dictstructures[str(mat)].append(mat)
            else:
                pass

        if save_to_file and (not fetch_results_only) and (path_to_save is not None):
            outputdict = {}
            outputdict['list_mat_id']=list_mat_id
            outputdict['structure_matching'] = dictstructures
            outputdict['additional_info'] = information_mat
            self._save_results_to_file(outputdict, path_to_save)

        if fetch_results_only:
            outputdict=self._get_precomputed_results(path_to_save)
            if not set(outputdict['list_mat_id'])==set(list_mat_id):
                #TODO: maybe different error?
                raise ValueError



            #checke, ob liste gleich ist!

        return outputdict

    def _print_to_screen_similar_structures(self,dict_to_print):
        for key,items in OrderedDict(sorted(dict_to_print['structure_matching'].items(), key=lambda t: len(t[1]), reverse=True)).items():
            if len(items)==1:
                print(key+' ('+str(dict_to_print['additional_info'][key])+', '+str(len(items)) +' representative): ', end='')
            else:
                print(key+' ('+str(dict_to_print['additional_info'][key])+', '+str(len(items)) +' representatives): ', end='')
            for item in items:
                print(item+' ('+str(dict_to_print['additional_info'][item])+'), ' ,end='')
            print()




# most frequent environment: this is also something that needs to be computed
# might need to be included in PaulingRules

class Pauling1OverAllAnalysis(OverAllAnalysis):


    def run(self, start_from_results=False, save_result_data=True,save_structure_analysis=True, restart_from_saved_structure_analyisis=False,path_to_save='Results/Results_First_Rule.json'):

        if not start_from_results:
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
            # TODO: save the analysis correctly
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source, save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split(
                                                                                 '.')[0] + "_structural_exceptions.json",fetch_results_only=restart_from_saved_structure_analyisis)
            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source, save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split(
                                                                                 '.')[0] + "_structurres_fulfilling.json",fetch_results_only=restart_from_saved_structure_analyisis)

            print("Exceptions")
            self._print_to_screen_similar_structures(dict_similarstructures_exceptions)
            print("Structures fulfilling the rule:")
            self._print_to_screen_similar_structures(dict_similarstructures_fulfilling)

        #print('Structures that fulfill rule')
        #print(self.structures_fulfillingrule)
        #print('Structures that do not fulfill rule')
        #print(self.structures_exceptions)
        #print('Structures that have environments/cations not covered')
        #print(self.structures_cannot_be_evaluated)
        #print(self.Plot_PSE_DICT)

    def _new_setup(self):
        # could include: newsetup or from file here!
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        self.Plot_PSE_DICT = {}
        # valence dependency can be introduced later

        for mat in list_mat[0:100]:
            print(mat)
            lse = self._get_lse_from_folder(mat, source=self.source)

            pauling1 = Pauling1(lse=lse, filenameradii='../univalent_cat_radii.json', onlylowerlimit=False)

            # one could also put this in a method!
            try:
                if pauling1.is_fulfilled():
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)
                print('Cannot be analyzed')

            # get details and add them
            Details = pauling1.get_details()['cat_dependency']
            self._add_dict_cat_dependency(self.Plot_PSE_DICT, Details)


class Pauling2OverAllAnalysis(OverAllAnalysis):

    def run(self, start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,save_structure_analysis=True, path_to_save='Results/Results_Second_Rule.json'):

        if not start_from_results:
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

        self._secondrule_plot(arraydev=self.arraydev_share * 100.0, relativefreqarray=self.relativefrequency * 100.0,
                              tot_stddev=self.tot_stddev * 100.0)

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
            # TODO: think how to load it from file in a save way
            dict_similarstructures_exceptions = self._get_similar_structures(self.structures_exceptions,
                                                                             source=self.source, save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structural_exceptions.json",fetch_results_only=restart_from_saved_structure_analyisis)

            dict_similarstructures_fulfilling = self._get_similar_structures(self.structures_fulfillingrule,
                                                                             source=self.source, save_to_file=save_structure_analysis,
                                                                             path_to_save=path_to_save.split('.')[
                                                                                              0] + "_structurres_fulfilling.json",fetch_results_only=restart_from_saved_structure_analyisis)

            #TODO: write a file with the output
            print('Exceptions:')
            self._print_to_screen_similar_structures(dict_similarstructures_exceptions)
            print('Structures fulfilling the rule')
            self._print_to_screen_similar_structures(dict_similarstructures_fulfilling)

    def stddev(self, lst, mean):
        """percental standard deviation of a sample"""
        sum = 0.0
        mn = mean
        for i in range(len(lst)):
            sum += pow((lst[i] - mn), 2)
        return np.sqrt([sum / (len(lst) - 1)]) / mn

    def stddev_normal(self, lst, mean):
        """?"""
        """standard deviation of a sample"""
        sum = 0.0
        mn = mean
        for i in range(len(lst)):
            sum += pow((lst[i] - mn), 2)
        return np.sqrt([sum / (len(lst) - 1)])

    def _secondrule_plot(self, arraydev, relativefreqarray, tot_stddev, maxpercentage=70):
        # TODO: include save to file
        font = {'size': 10}

        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        matplotlib.rc('font', **font)

        font = {'family': 'normal',
                'size': 22}

        matplotlib.rc('font', **font)
        plt.plot((arraydev), relativefreqarray)
        plt.plot([tot_stddev, tot_stddev], [0, 100.0], linewidth=3.0)
        plt.ylim([0, 100.0])
        plt.xlim([0, maxpercentage])
        plt.xticks(range(0, maxpercentage + 10, 10))
        plt.yticks(range(0, 100 + 20, 20))
        plt.xlabel("Absolute Deviation (%) from the ideal valence -2")
        plt.ylabel("Oxygen Atoms (%)")
        # plt.savefig('peroxygen.svg')
        plt.show()

    def _new_setup(self):
        list_mat = self._get_list_materials(source=self.source, onlybinaries=self.onlybinaries)

        self.structures_fulfillingrule = []
        self.structures_exceptions = []
        self.structures_cannot_be_evaluated = []

        array_bvs = []
        self.Plot_PSE_DICT = {}
        self.present_env = {}
        # valence dependency can be introduced later

        for mat in list_mat:
            lse = self._get_lse_from_folder(mat, source=self.source)

            pauling2 = Pauling2(lse=lse)

            # one could also put this in a method!
            try:
                if pauling2.is_fulfilled():
                    self.structures_fulfillingrule.append(mat)
                else:
                    self.structures_exceptions.append(mat)
            except RuleCannotBeAnalyzedError:
                self.structures_cannot_be_evaluated.append(mat)

            # get details and add them
            Details = pauling2.get_details()
            bvs = Details['bvs_for_each_anion']

            self._add_dict_cat_dependency(self.Plot_PSE_DICT, Details['elementwise_fulfillment'])

            # add cations to each other, to count total number of cations
            self._add_dict_cat_dependency(self.present_env, Details['cations_in_structure'],
                                          number_of_elements_to_add=1)

            array_bvs.extend(bvs)

            # print(self.array_bvs)

        # alles stueck fuer stueck durchschauen - vor allem den plot!
        self.bs_sum_mean = np.mean(array_bvs)

        self.tot_stddev = self.stddev(array_bvs, np.mean(array_bvs))
        ideal_bs = 2.0
        bs_dev = [abs(x - ideal_bs) for x in array_bvs]
        # print(bs_dev)
        # print(bs_sum_mean)
        # print(tot_stddev)
        self.frequency = []
        dev_array = np.arange(0, max(bs_dev), 0.01)
        for step in dev_array:
            self.frequency.append(len([x for x in bs_dev if x < step]))

        # print(self.frequency)
        # print(self.dev_array)

        self.arraydev_share = np.array(dev_array) / ideal_bs
        self.relativefrequency = np.array(self.frequency) / float(len(array_bvs))

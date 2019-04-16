from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, \
    Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, Pauling1Entropy, Pauling1Frequency, AllPaulingOverAllAnalysis

from Classes_for_statistics import OverAllAnalysis

#
# over=OverAllAnalysis()
# print(over._get_similar_structures(["mp-7000","mp-5986", "mp-306"],source='experimental',save_to_file=False,fetch_results_only=False,start_from_Matching=True,restart_from_matching=False))


# TODO: write new allmaterials list! -> has some problems -> some duplicates, which are problematic
# TODO: clean this
#
#
# newclass=Pauling1OverAllAnalysis(source='experimental',onlybinaries=False,plot_element_dependend_analysis=True,analyse_structures=False,use_prematching=True,lowest_number_environments_for_plot=50)
# newclass.run(start_from_results=False,save_result_data=True,save_structure_analysis=True,restart_from_saved_structure_analyisis=False, path_to_save="Results/Results_First_Rule_ICSD.json")

# newclass=Pauling1Frequency(source='experimental',onlybinaries=False,plot_element_dependend_analysis=True, lowest_number_environments_for_plot=50,lower_limit_plot=0.0, upper_limit_plot=1.0)
# newclass.run(start_from_results=False,save_result_data=True,path_to_save='Results/Results_First_ICSD_Rule_Most_Frequent.json')
#
#
# newclass=Pauling1Entropy(source='experimental',onlybinaries=False,plot_element_dependend_analysis=True, lowest_number_environments_for_plot=50,lower_limit_plot=0.1, upper_limit_plot=1.0)
# newclass.run(start_from_results=False,save_result_data=True,path_to_save='Results/Results_First_ICSD_Rule_Entropy.json')


# newclass = Pauling2OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,lowest_number_environments_for_plot=50, lower_limit_plot=0.1,upper_limit_plot=0.8,analyse_structures=False,use_prematching=True)
# newclass.run(start_from_results=False, save_result_data=True,path_to_save='Results/Results_Second_Rule_ICSD.json',save_structure_analysis=True,restart_from_saved_structure_analyisis=False)


newclass = Pauling3OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections_exp',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule_ICSD.json')
#
# newclass = Pauling4OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.8, upper_limit_plot=1.0,
#                                    analyse_structures=True, use_prematching=True)
# newclass.run(show_plot=True, start_from_connections=True, save_connections=True,
#              connections_folder='AnalysisConnections_exp',
#              start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
#              save_structure_analysis=True, path_to_save='Results/Results_Fourth_Rule_ICSD.json')


# newclass = Pauling5OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
#                                    analyse_structures=True, use_prematching=True)
#
# newclass.run(show_plot=True, remove_elements_low_entropy=False, threshold_remove_elements=0.9,
#              start_from_connections=True, save_connections=True,
#              connections_folder='AnalysisConnections_5thRule_exp',
#              start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
#              save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule_ICSD.json')

# remove elements with low entropy!
# newclass = Pauling5OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
#                                    analyse_structures=True, use_prematching=True)
#
# newclass.run(show_plot=True, remove_elements_low_entropy=True, threshold_remove_elements=0.9,
#              start_from_connections=True, save_connections=True,
#              connections_folder='AnalysisConnections_5thRule_exp',
#              start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
#              save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule_ICSD_remove_low_entropy.json')
#
# newclass = AllPaulingOverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
#                                      lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
#                                      analyse_structures=True, use_prematching=True)
#
# newclass.run(remove_elements_low_entropy=False, start_from_connections=True, save_connections=True,
#              connections_folder34='AnalysisConnections_exp', connections_folder5='AnalysisConnections_5thRule_exp',
#              start_from_results=False, save_result_data=True,
#              restart_from_saved_structure_analyisis=False, save_structure_analysis=True,
#              path_to_save='Results/Results_AllRules_ICSD.json', start_material=None, stop_material=None,
#              threshold_remove_elements=0.95)



# #To Match all structures:
# MatchAllStructures(source='experimental',startmaterial=0,stopmaterial=10,restart_from_matching=False)
# MatchAllStructures(source='experimental',startmaterial=10,stopmaterial=20,restart_from_matching=True)

# MatchAllStructures(source='experimental',startmaterial=3000,stopmaterial=5891,restart_from_matching=True)


# TODO: make different result folders: normal results, binaries, only very symmetric, ICSD

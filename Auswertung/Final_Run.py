from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, \
    Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, MatchAllStructures, Pauling1Entropy, Pauling1Frequency

from Classes_for_statistics import OverAllAnalysis

#
# over=OverAllAnalysis()
# print(over._get_similar_structures(["mp-7000","mp-5986", "mp-306"],source='MP',save_to_file=False,fetch_results_only=False,start_from_Matching=True,restart_from_matching=False))


# TODO: write new allmaterials list! -> has some problems -> some duplicates, which are problematic
# TODO: clean this
#
#
# newclass=Pauling1OverAllAnalysis(source='MP',onlybinaries=False,plot_element_dependend_analysis=True,analyse_structures=True,use_prematching=True)
# newclass.run(start_from_results=False,save_result_data=True,save_structure_analysis=True,restart_from_saved_structure_analyisis=False)

# newclass=Pauling1Frequency(source='MP',onlybinaries=False,plot_element_dependend_analysis=True, lowest_number_environments_for_plot=50,lower_limit_plot=0.0, upper_limit_plot=1.0)
# newclass.run(start_from_results=False,save_result_data=True,path_to_save='Results/Results_First_Rule_Most_Frequent.json')


# newclass=Pauling1Entropy(source='MP',onlybinaries=False,plot_element_dependend_analysis=True, lowest_number_environments_for_plot=50,lower_limit_plot=0.55, upper_limit_plot=1.0)
# newclass.run(start_from_results=False,save_result_data=True,path_to_save='Results/Results_First_Rule_Entropy.json')


# newclass = Pauling2OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,lowest_number_environments_for_plot=50, lower_limit_plot=0.1,upper_limit_plot=0.7,analyse_structures=True,use_prematching=True)
# newclass.run(start_from_results=False, save_result_data=True,path_to_save='Results/Results_Second_Rule.json',save_structure_analysis=True,restart_from_saved_structure_analyisis=False)
#
# # TODO: third rule and fourth rule!
# # TODO: make it possible to restart fourth rule from connections from third rule
# # save_connection_data as well as save_result_data
# # start_from_connections as well as start_from_results
# # has to run very long on a computer in the computing center
# newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
#                                    analyse_structures=True,use_prematching=True)
# newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
#              start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
#              save_structure_analysis=True, path_to_save='Results/Results_Third_Rule.json')
#
# newclass = Pauling4OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=False,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
#                                    analyse_structures=True,use_prematching=True)
# newclass.run(show_plot=False,start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
#              start_from_results=True, save_result_data=True, restart_from_saved_structure_analyisis=True,
#              save_structure_analysis=True, path_to_save='Results/Results_Fourth_Rule.json')
#
#

newclass = Pauling5OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=False, use_prematching=True)

newclass.run(show_plot=True,remove_elements_low_entropy=True,threshold_remove_elements=0.90, start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule.json')

# #To Match all structures:
# MatchAllStructures(source='MP',startmaterial=0,stopmaterial=10,restart_from_matching=False)
# MatchAllStructures(source='MP',startmaterial=10,stopmaterial=20,restart_from_matching=True)

# MatchAllStructures(source='MP',startmaterial=3000,stopmaterial=5891,restart_from_matching=True)


# TODO: make different result folders: normal results, binaries, only very symmetric, ICSD

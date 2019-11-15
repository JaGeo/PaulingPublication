from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, \
    Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, Pauling1Entropy, AllPaulingOverAllAnalysis, \
    AllPaulingOverAllAnalysis_Final_Summary, HowMany, EntropyDeviationFrom2ndRuleDiagram

print("Calculates how many environments of a certain element are present in the data set")
newclass = HowMany(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                   lower_limit_plot=0., upper_limit_plot=2000.0, lowest_number_of_environments_considered=0,
                   upper_number_of_environments_considered=2000.0)
newclass.run(start_from_results=False, save_result_data=True, path_to_save='Results/Results_Count_Compounds.json')

print("Calculates how many environments of a certain element are present in the data set")
newclass = HowMany(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                   lower_limit_plot=2000., upper_limit_plot=8500.0, lowest_number_of_environments_considered=2000,
                   upper_number_of_environments_considered=8500)
newclass.run(start_from_results=True, save_result_data=True, path_to_save='Results/Results_Count_Compounds.json')

print("Analysis of the first rule")
print("Check with the help of the univalent radii")
newclass = Pauling1OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   analyse_structures=True, use_prematching=True,
                                   lowest_number_environments_for_plot=50)
newclass.run(start_from_results=False, save_result_data=True, save_structure_analysis=True,
             restart_from_saved_structure_analysis=False)

print("Evaluation of Shannon entropy for the coordination environments of each element")
newclass = Pauling1Entropy(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                           lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=1.0)
newclass.run(start_from_results=False, save_result_data=True, path_to_save='Results/Results_First_Rule_Entropy.json')

# newclass= EntropyDeviationFrom2ndRuleDiagram(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
#                             lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=1.0)
# newclass.run(start_from_results=False, save_result_data=True, path_to_save='Results/Results_Second_Rule_Entropy_vs_Deviation.json')

print("Analysis of the 2nd rule: all coordination environments")
newclass = Pauling2OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=0.8,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_results=False, save_result_data=True, path_to_save='Results/Results_Second_Rule.json',
             save_structure_analysis=True, restart_from_saved_structure_analysis=False, show_histogram=True,
             stepsize_histogram=0.1)

print(
    "Analysis of the 2nd rule: all coordination environments; not only the environments with the highest probability are considered but all of them until a fraction of 0.3")
newclass = Pauling2OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=0.8,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_Second_Rule_optimized_environments.json',
             save_structure_analysis=True, restart_from_saved_structure_analysis=False, show_histogram=True,
             stepsize_histogram=0.1, optimized_environments=True)


print("Analysis of the 2nd rule: only very symmetric coordination environments")
newclass = Pauling2OverAllAnalysis(source='MP_very_symmetric', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.5, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_Second_Rule_very_symmetric.json', save_structure_analysis=True,
             restart_from_saved_structure_analysis=False)

print("Analysis of third rule: all coordination numbers")
newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule.json')

print(
    "Analysis of third rule: all coordination numbers. In the elementwise analysis, also edges will be treated as exceptions.")
newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.45, upper_limit_plot=0.85,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=False,
             path_to_save='Results/Results_Third_Rule_edges_as_exceptions_in_elementwise_analysis.json',
             EdgesAsAdditionalExceptions=True)

print("Analysis of third rule: only coordination numbers smaller or equal to 8")
newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.86, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule_max_CN.json', maxCN=8)

print(
    "Analysis of the third rule: only structures with ALL CN smaller or equal to 8 will be considered! This is different from the last analysis!")
newclass = AllPaulingOverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                     lowest_number_environments_for_plot=25, lower_limit_plot=0.0, upper_limit_plot=0.3,
                                     analyse_structures=True, use_prematching=True)
newclass.run(remove_elements_low_entropy=False, start_from_connections=True,
             save_connections=True, connections_folder34='AnalysisConnections',
             connections_folder5='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True,
             restart_from_saved_structure_analysis=False, save_structure_analysis=True,
             path_to_save='Results/Results_Third_Rule_max_CN_for_whole_structure.json', threshold_remove_elements=0.90,
             start_material=None,
             stop_material=None, adapt_first_fourth_and_fifth_rules=True, ignore_first_rule=True,
             ignore_second_rule=True, ignore_third_rule=False, ignore_fourth_rule=True, ignore_fifth_rule=True,
             remove_structures_with_CN_larger_8=True)


print("Analysis of fourth rule")
newclass = Pauling4OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.8, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(show_plot=True, start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fourth_Rule.json')

print("Analysis of fifth rule: no Shannon entropy considered")
newclass = Pauling5OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)

newclass.run(show_plot=True, remove_elements_low_entropy=False, threshold_remove_elements=0.9,
             start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule.json')

# #remove elements with low entropy!
print("Analysis of fifth rule: remove cations with low Shannon entropy from analysis")
newclass = Pauling5OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)

newclass.run(show_plot=True, remove_elements_low_entropy=True, threshold_remove_elements=0.9,
             start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule_remove_low_entropy.json')

print("Analyse Pauling rules 2-5 and print structural exceptions")
newclass = AllPaulingOverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
                                     lowest_number_environments_for_plot=25, lower_limit_plot=0.0, upper_limit_plot=0.3,
                                     analyse_structures=True, use_prematching=True)
newclass.run(remove_elements_low_entropy=False, start_from_connections=True,
             save_connections=True, connections_folder34='AnalysisConnections',
             connections_folder5='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True,
             restart_from_saved_structure_analysis=False, save_structure_analysis=True,
             path_to_save='Results/Results_AllRules.json', threshold_remove_elements=0.90, start_material=None,
             stop_material=None, adapt_first_fourth_and_fifth_rules=True, ignore_first_rule=True,
             ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False, ignore_fifth_rule=False,
             remove_structures_with_CN_larger_8=False)

print(
    "Analyse Pauling rules 2-5 and test the influence of each of the rules. Adapt the criteria for structures to assess the fourth and fifth rule.")

newclass = AllPaulingOverAllAnalysis_Final_Summary(source='MP', onlybinaries=False,
                                                   plot_element_dependend_analysis=False,
                                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6,
                                                   upper_limit_plot=1.0,
                                                   analyse_structures=False, use_prematching=True)

newclass.run(remove_elements_low_entropy=False, start_from_connections=True, save_connections=True,
             connections_folder34='AnalysisConnections', connections_folder5='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_AllRules_final_plot.json', start_material=None, stop_material=None,
             threshold_remove_elements=0.90)

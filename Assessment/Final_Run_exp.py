from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, \
    Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, Pauling1Entropy, Pauling1Frequency, AllPaulingOverAllAnalysis, \
    AllPaulingOverAllAnalysis_Final_Summary

print("Analysis of the first rule")
print("Check with the help of the univalent radii")
newclass = Pauling1OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   analyse_structures=False, use_prematching=True,
                                   lowest_number_environments_for_plot=50)
newclass.run(start_from_results=False, save_result_data=True, save_structure_analysis=True,
             restart_from_saved_structure_analysis=False, path_to_save="Results/Results_First_Rule_exp.json")

print("Evaluation of Shannon entropy for the coordination environments of each element")
newclass = Pauling1Entropy(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                           lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=1.0)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_First_exp_Rule_Entropy.json')

print("Analysis of the 2nd rule")

newclass = Pauling2OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=0.8,
                                   analyse_structures=False, use_prematching=True)
newclass.run(start_from_results=False, save_result_data=True, path_to_save='Results/Results_Second_Rule_exp.json',
             save_structure_analysis=True, restart_from_saved_structure_analysis=False)

print("Analysis of third rule: all coordination numbers")
newclass = Pauling3OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=False, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections_exp',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule_exp.json')

print("Analysis of third rule: only coordination numbers smaller or equal to 8")
newclass = Pauling3OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.83, upper_limit_plot=1.0,
                                   analyse_structures=False, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections_exp',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule_max_CN_exp.json', maxCN=8)

print("Analysis of fourth rule")
newclass = Pauling4OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.7, upper_limit_plot=1.0,
                                   analyse_structures=False, use_prematching=True)
newclass.run(show_plot=True, start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_exp',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fourth_Rule_exp.json')

print("Analysis of fifth rule: no Shannon entropy considered")
newclass = Pauling5OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=False, use_prematching=True)

newclass.run(show_plot=True, remove_elements_low_entropy=False, threshold_remove_elements=0.9,
             start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule_exp',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule_exp.json')
#
print("Analysis of fifth rule: remove cations with low Shannon entropy from analysis")
# remove elements with low entropy!
newclass = Pauling5OverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=False, use_prematching=True)

newclass.run(show_plot=True, remove_elements_low_entropy=True, threshold_remove_elements=0.9,
             start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule_exp',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule_exp_remove_low_entropy.json')

print("Analyse Pauling rules 2-5 and print percentage of structures fulfilling the rules")
newclass = AllPaulingOverAllAnalysis(source='experimental', onlybinaries=False, plot_element_dependend_analysis=True,
                                     lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=0.3,
                                     analyse_structures=False, use_prematching=True)
newclass.run(remove_elements_low_entropy=False, start_from_connections=True,
             save_connections=True, connections_folder34='AnalysisConnections_exp',
             connections_folder5='AnalysisConnections_5thRule_exp',
             start_from_results=False, save_result_data=True,
             restart_from_saved_structure_analysis=False, save_structure_analysis=True,
             path_to_save='Results/Results_AllRules_exp.json', threshold_remove_elements=0.90, start_material=None,
             stop_material=None, adapt_first_fourth_and_fifth_rules=True, ignore_first_rule=True,
             ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False, ignore_fifth_rule=False,
             remove_structures_with_CN_larger_8=False)

print(
    "Analyse Pauling rules 2-5 and test the influence of each of the rules. Adapt the criteria for structures to assess the fourth and fifth rule.")

newclass = AllPaulingOverAllAnalysis_Final_Summary(source='experimental', onlybinaries=False,
                                                   plot_element_dependend_analysis=False,
                                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6,
                                                   upper_limit_plot=1.0,
                                                   analyse_structures=False, use_prematching=True)

newclass.run(remove_elements_low_entropy=False, start_from_connections=True, save_connections=True,
             connections_folder34='AnalysisConnections_exp', connections_folder5='AnalysisConnections_5thRule_exp',
             start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_AllRules_final_plot_exp.json', start_material=None, stop_material=None,
             threshold_remove_elements=0.90)


from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, \
    Pauling4OverAllAnalysis, Pauling5OverAllAnalysis, Pauling1Entropy, AllPaulingOverAllAnalysis, \
    AllPaulingOverAllAnalysis_Final_Summary, HowMany

print("Only binaries are analysed in the following:")
print("Calculates how many environments of a certain element are present in the data set")
newclass = HowMany(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                   lower_limit_plot=0., upper_limit_plot=150.0, lowest_number_of_environments_considered=0,
                   upper_number_of_environments_considered=150.0)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_Count_Compounds_binaries.json')

newclass = HowMany(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                   lower_limit_plot=1250., upper_limit_plot=1500.0, lowest_number_of_environments_considered=1250,
                   upper_number_of_environments_considered=1500.0)
newclass.run(start_from_results=True, save_result_data=True,
             path_to_save='Results/Results_Count_Compounds_binaries.json')

print("Analysis of the first rule")
print("Check with the help of the univalent radii")
newclass = Pauling1OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   analyse_structures=True, use_prematching=True,
                                   lowest_number_environments_for_plot=50)
newclass.run(start_from_results=False, save_result_data=True, save_structure_analysis=True,
             restart_from_saved_structure_analysis=False, path_to_save='Results/Results_First_Rule_only_binaries.json')

print("Evaluation of Shannon entropy for the coordination environments of each element")
newclass = Pauling1Entropy(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                           lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=1.0)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_First_Rule_Entropy_only_binaries.json')

print("Analysis of the 2nd rule: all environments")
newclass = Pauling2OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.1, upper_limit_plot=0.8,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_Second_Rule_only_binaries.json', save_structure_analysis=True,
             restart_from_saved_structure_analysis=False)

print("Analysis of the 2nd rule: only very symmetric coordination environments")
newclass = Pauling2OverAllAnalysis(source='MP_very_symmetric', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.5, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_Second_Rule_very_symmetric_only_binaries.json', save_structure_analysis=True,
             restart_from_saved_structure_analysis=False)

print("Analysis of third rule: all coordination numbers")
newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule_only_binaries.json')

print("Analysis of third rule: only coordination numbers smaller or equal to 8")
newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule_max_CN_only_binaries.json', maxCN=8)

print("Analysis of fourth rule")
newclass = Pauling4OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.8, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)
newclass.run(show_plot=True, start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fourth_Rule_only_binaries.json')

print("Analysis of fifth rule: no Shannon entropy considered")
newclass = Pauling5OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)

newclass.run(show_plot=True, remove_elements_low_entropy=False, threshold_remove_elements=0.9,
             start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule_only_binaries.json')

print("Analysis of fifth rule: remove cations with low Shannon entropy from analysis")
newclass = Pauling5OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                   analyse_structures=True, use_prematching=True)

newclass.run(show_plot=True, remove_elements_low_entropy=True, threshold_remove_elements=0.9,
             start_from_connections=True, save_connections=True,
             connections_folder='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analysis=False,
             save_structure_analysis=True,
             path_to_save='Results/Results_Fifth_Rule_remove_low_entropy_only_binaries.json')

print("Analyse Pauling rules 2-5 and print structural exceptions")
newclass = AllPaulingOverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                     lowest_number_environments_for_plot=50, lower_limit_plot=0.6, upper_limit_plot=1.0,
                                     analyse_structures=True, use_prematching=True)
newclass.run(remove_elements_low_entropy=False, start_from_connections=True,
             save_connections=True, connections_folder34='AnalysisConnections',
             connections_folder5='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True,
             restart_from_saved_structure_analysis=False, save_structure_analysis=True,
             path_to_save='Results/Results_AllRules_binaries.json', threshold_remove_elements=0.90, start_material=None,
             stop_material=None, adapt_first_fourth_and_fifth_rules=True, ignore_first_rule=True,
             ignore_second_rule=False, ignore_third_rule=False, ignore_fourth_rule=False, ignore_fifth_rule=False,
             remove_structures_with_CN_larger_8=False)

print(
    "Analyse Pauling rules 2-5 and test the influence of each of the rules. Adapt the criteria for structures to assess the fourth and fifth rule.")

newclass = AllPaulingOverAllAnalysis_Final_Summary(source='MP', onlybinaries=True,
                                                   plot_element_dependend_analysis=False,
                                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.6,
                                                   upper_limit_plot=1.0,
                                                   analyse_structures=False, use_prematching=True)

newclass.run(remove_elements_low_entropy=False, start_from_connections=True, save_connections=True,
             connections_folder34='AnalysisConnections', connections_folder5='AnalysisConnections_5thRule',
             start_from_results=False, save_result_data=True,
             path_to_save='Results/Results_AllRules_final_plot_binaries.json', start_material=None, stop_material=None,
             threshold_remove_elements=0.90)

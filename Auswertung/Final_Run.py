from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis

#newclass=Pauling1OverAllAnalysis(source='MP',onlybinaries=False,plot_element_dependend_analysis=True)
#newclass.run(start_from_results=False,save_result_data=False,save_structure_analysis=True,restart_from_saved_structure_analyisis=True)


newclass = Pauling2OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,lowest_number_environments_for_plot=10, lower_limit_plot=0.0,upper_limit_plot=1.0,analyse_structures=False)
newclass.run(start_from_results=False, save_result_data=False,path_to_save='Results/Results_Second_Rule.json',save_structure_analysis=False,restart_from_saved_structure_analyisis=False)

#TODO: third rule and fourth rule!
#TODO: make it possible to restart fourth rule from connections from third rule

#TODO: make different result folders: normal results, binaries, only very symmetric, ICSD
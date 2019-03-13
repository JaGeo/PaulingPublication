from Classes_for_statistics import Pauling1OverAllAnalysis, Pauling2OverAllAnalysis, Pauling3OverAllAnalysis, Pauling4OverAllAnalysis, Pauling5OverAllAnalysis

# newclass=Pauling1OverAllAnalysis(source='MP',onlybinaries=False,plot_element_dependend_analysis=True)
# newclass.run(start_from_results=False,save_result_data=False,save_structure_analysis=True,restart_from_saved_structure_analyisis=True)


# newclass = Pauling2OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,lowest_number_environments_for_plot=10, lower_limit_plot=0.0,upper_limit_plot=1.0,analyse_structures=False)
# newclass.run(start_from_results=False, save_result_data=False,path_to_save='Results/Results_Second_Rule.json',save_structure_analysis=False,restart_from_saved_structure_analyisis=False)

# TODO: third rule and fourth rule!
# TODO: make it possible to restart fourth rule from connections from third rule
# save_connection_data as well as save_result_data
# start_from_connections as well as start_from_results
# has to run very long on a computer in the computing center
newclass = Pauling3OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
                                   lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
                                   analyse_structures=True)
newclass.run(start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
             start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
             save_structure_analysis=True, path_to_save='Results/Results_Third_Rule.json') #, start_material=0,
             #stop_material=19)


#maybe run the binaries during the night
# newclass = Pauling4OverAllAnalysis(source='MP', onlybinaries=True, plot_element_dependend_analysis=True,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
#                                    analyse_structures=True)
# newclass.run(show_plot=False,start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections',
#              start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
#              save_structure_analysis=True, path_to_save='Results/Results_Fourth_Rule.json')


# newclass = Pauling5OverAllAnalysis(source='MP', onlybinaries=False, plot_element_dependend_analysis=True,
#                                    lowest_number_environments_for_plot=50, lower_limit_plot=0.0, upper_limit_plot=1.0,
#                                    analyse_structures=True)
# newclass.run(show_plot=False,start_from_connections=True, save_connections=True, connections_folder='AnalysisConnections_5thRule',
#              start_from_results=False, save_result_data=True, restart_from_saved_structure_analyisis=False,
#              save_structure_analysis=True, path_to_save='Results/Results_Fifth_Rule.json',start_material=None,stop_material=20)




# TODO: make different result folders: normal results, binaries, only very symmetric, ICSD

from PaulingRules import Pauling2_optimized_environments, Pauling2
from Classes_for_statistics import OverAllAnalysis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import json
# def test_recursive(list1,list2):
#     new_list=[]
#     for el in list1:
#         for el2 in list2:
#             new=el.copy()
#             new.append(el2)
#             new_list.append(new)
#     return new_list
#
#
# whole_list=[[0],[0,1],[0,1,2]]
#
# output=[whole_list[0]]
# for part_list in whole_list[1:]:
#     output=test_recursive(output,part_list)
#

# print(test_recursive(list1,list2))
#
# # #TODO: search for a well-suited example
analysis = OverAllAnalysis()
#set of materials that shows problems in analysis of second rule
list_lse = analysis._get_list_materials()
# new_lse_list=[]
# for mat in list_lse:
#     if mat not in ['mp-559234', 'mp-558239', 'mp-559460', 'mp-555131', 'mp-565811', 'mp-31358', 'mp-558128', 'mp-32131', 'mp-652260', 'mp-561433', 'mp-18136', 'mp-21099', 'mp-560516', 'mp-19035', 'mp-510281', 'mp-656144', 'mp-19560', 'mp-541454', 'mp-680384', 'mp-566937', 'mp-579382', 'mp-17605', 'mp-550579', 'mp-18977', 'mp-6760', 'mp-645568', 'mp-554929', 'mp-19093', 'mp-555387', 'mp-6702', 'mp-19442', 'mp-561689', 'mp-10417', 'mp-23352', 'mp-28708', 'mp-4112', 'mp-18992', 'mp-560605', 'mp-541436', 'mp-505110', 'mp-543034', 'mp-561634', 'mp-15379', 'mp-19400', 'mp-561055', 'mp-565692', 'mp-31624', 'mp-555515', 'mp-24854', 'mp-9480', 'mp-566311', 'mp-568561', 'mp-17770', 'mp-19470', 'mp-566249', 'mp-653182', 'mp-25053', 'mp-25152', 'mp-6318', 'mp-545735', 'mp-554603', 'mp-559373', 'mp-4396', 'mp-3919', 'mp-631631', 'mp-19586', 'mp-510624', 'mp-18986', 'mp-557382', 'mp-554962', 'mp-557997', 'mp-18750', 'mp-647385', 'mp-616597', 'mp-557372']:
#         new_lse_list.append(mat)
#
# with open("allmaterials.json",'w') as f:
#     json.dump({"is_clear_materials":new_lse_list},f)

#     #
#     lse = analysis._get_lse_from_folder(mat)
#     # print(lse.structure)
#     # print(lse.structure)
#     spanalyzer = SpacegroupAnalyzer(structure=lse.structure)
#     equivalent_indices = spanalyzer.get_symmetrized_structure().equivalent_indices
#     print(equivalent_indices)
#
#     lgf = LocalGeometryFinder()
#
#     lgf.setup_structure(structure=lse.structure)
#     se = lgf.compute_structure_environments(only_cations=True,valences=lse.valences)
#
#     from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
#     strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
#
#     lse2 = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
#
#
#     #for ice,ce in enumerate(lse.coordination_environments):
#     for equi in equivalent_indices:
#         print("Set1")
#         for i in equi:
#             if lse.coordination_environments[i]!=None:
#                 print(lse.coordination_environments[i][0]['ce_symbol'])
#         print("Set2")
#         for i in equi:
#             if lse2.coordination_environments[i] != None:
#                 print(lse2.coordination_environments[i][0]['ce_symbol'])
#
#
# exit()
#
# lgf = LocalGeometryFinder()
#
# lgf.setup_structure(structure=lse.structure)
# se = lgf.compute_structure_environments(only_cations=True,valences=lse.valences)
#
# from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
# strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
#
# lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)
# for ce in lse.coordination_environments:
#     if ce is not None:
#         print(len(ce))
#         print(ce)

# pauling2=Pauling2_optimized_environments(lse=lse,perc=0.1)
# print(pauling2.is_fulfilled())
# exit()
#
list_lse = analysis._get_list_materials()

fulfilled = 0
all = 0
exceptions = []
for mat in list_lse:
    try:
    # if not mat in ["mp-17304","mp-15472","mp-581229","mp-680384","mp-706588","mp-581345","mp-680432","mp-561374",'mp-557970', 'mp-5637', 'mp-645568', 'mp-19128', 'mp-554219', 'mp-556195', 'mp-653973', 'mp-20546', 'mp-565522', 'mp-12010', 'mp-19660', 'mp-560919','mp-557970', 'mp-5637', 'mp-645568', 'mp-19128', 'mp-554219','mp-566537', 'mp-10444',"mp-581644", 'mp-556195', 'mp-653973', 'mp-20546', 'mp-565522', 'mp-12010', 'mp-19660', 'mp-560919',"mp-3974"]:
        print(mat)
        # try:
        lse = analysis._get_lse_from_folder(mat=mat, source="MP")
        pauling2 = Pauling2_optimized_environments(lse=lse, perc=0.3)
        if (pauling2.is_fulfilled()):
            fulfilled += 1
        all += 1
    except:
        exceptions.append(mat)

print(float(fulfilled) / float(all))
print(fulfilled)
print(all)
print(len(list_lse))
print(exceptions)
# TODO: check if there is a dependency on the materials list order somewhere

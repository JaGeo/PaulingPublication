from PaulingRules import Pauling2_optimized_environments, Pauling2
from Classes_for_statistics import OverAllAnalysis
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments

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

# #TODO: search for a well-suited example
analysis = OverAllAnalysis()
for mat in ['mp-17304', 'mp-557970', 'mp-20546', 'mp-15472']:
    #
    lse = analysis._get_lse_from_folder(mat)
    # print(lse.structure)
    # print(lse.structure)
    spanalyzer = SpacegroupAnalyzer(structure=lse.structure)
    equivalent_indices = spanalyzer.get_symmetrized_structure().equivalent_indices
    print(equivalent_indices)

    lgf = LocalGeometryFinder()

    lgf.setup_structure(structure=lse.structure)
    se = lgf.compute_structure_environments(only_cations=True,valences=lse.valences)

    from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
    strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()

    lse2 = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)


    #for ice,ce in enumerate(lse.coordination_environments):
    for equi in equivalent_indices:
        print("Set1")
        for i in equi:
            if lse.coordination_environments[i]!=None:
                print(lse.coordination_environments[i][0]['ce_symbol'])
        print("Set2")
        for i in equi:
            if lse2.coordination_environments[i] != None:
                print(lse2.coordination_environments[i][0]['ce_symbol'])


exit()
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

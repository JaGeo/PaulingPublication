from PaulingRules import Pauling1, Pauling2, Pauling3, Pauling4, Pauling5, \
    RuleCannotBeAnalyzedError, is_an_oxide_and_no_env_for_O
from pymatgen.analysis.chemenv.coordination_environments.coordination_geometry_finder import LocalGeometryFinder
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.chemenv.coordination_environments.chemenv_strategies import MultiWeightsChemenvStrategy
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from pymatgen.analysis.bond_valence import BVAnalyzer

all_rules_can_be_evaluated = True

# first: get a structure, it has to be an oxide at the moment

a = MPRester()
#struct = a.get_structure_by_material_id('mp-1788')
struct = a.get_structure_by_material_id('mp-12236')

# get the valences, you can also insert them manually
bva = BVAnalyzer()
valences = bva.get_valences(structure=struct)

# Setup the local geometry finder
lgf = LocalGeometryFinder()
lgf.setup_structure(structure=struct)

# Get the StructureEnvironments
# only the environments of the cations should be computed, otherwise, the code will fail
se = lgf.compute_structure_environments(maximum_distance_factor=1.41, only_cations=True, valences=valences)

strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()

lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)

# test if you really have an oxide!
is_an_oxide_and_no_env_for_O(lse)

# # Test First Rule
# pauling1=Pauling1(lse=lse,onlylowerlimit=False)
#
# try:
#     pauling1_fulfilled=pauling1.is_fulfilled()
#     if pauling1_fulfilled:
#         print("Pauling's first rule is fulfilled")
#     else:
#         print("Pauling's first rule is not fulfilled")
# except RuleCannotBeAnalyzedError:
#     print("Pauling's first rule cannot be analyzed")
#     all_rules_can_be_evaluated=False
#
# # Test Second Rule
# pauling2=Pauling2(lse)
#
# try:
#     pauling2_fulfilled=pauling2.is_fulfilled()
#     if pauling2_fulfilled:
#         print("Pauling's second rule is fulfilled")
#     else:
#         print("Pauling's second rule is not fulfilled")
# except RuleCannotBeAnalyzedError:
#     print("Pauling's second rule cannot be analyzed")
#     all_rules_can_be_evaluated=False
#
# # Test Third Rule
# pauling3=Pauling3()
# pauling3.newsetup(lse,save_to_file=False)
#
#
# pauling3_fulfilled=pauling3.is_fulfilled()
# if pauling3_fulfilled:
#     print("Pauling's third rule is fulfilled")
# else:
#     print("Pauling's third rule is not fulfilled")
#
#
# # Test Fourth Rule
# pauling4=Pauling4()
# pauling4.newsetup(lse,save_to_file=False)
#
# try:
#     pauling4_fulfilled=pauling4.is_fulfilled()
#     if pauling4_fulfilled:
#         print("Pauling's fourth rule is fulfilled")
#     else:
#         print("Pauling's fourth rule is not fulfilled")
# except RuleCannotBeAnalyzedError:
#     print("Pauling's fourth rule cannot be analyzed")
#     all_rules_can_be_evaluated=False
#
#

# test Fifth Rule
leave_out_list = ["Er"]
pauling5 = Pauling5()
pauling5.newsetup(lse, save_to_file=False)
try:
    pauling5_fulfilled = pauling5.is_fulfilled(leave_out_list=leave_out_list)
    if pauling5_fulfilled:
        print("Pauling's fifth rule is fulfilled")
    else:
        print("Pauling's fifth rule is not fulfilled")
except RuleCannotBeAnalyzedError:
    print("Pauling's fifth rule cannot be analyzed")
    all_rules_can_be_evaluated = False

print(pauling5.get_details(options='CN', leave_out_list=leave_out_list))

#
# if all_rules_can_be_evaluated:
#     if pauling1_fulfilled and pauling2_fulfilled and pauling3_fulfilled and pauling4_fulfilled and pauling5_fulfilled:
#         print('All Rules are fulfilled.')

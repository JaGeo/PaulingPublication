from Graph_Test import Graph_Connection
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
from PaulingRules import Pauling1_general_limit_rule
import json


#test visualization


# # "mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-5986.json", "mp-306.json"
# with open("tests/mp-7000.json", "r") as f:
#     dict_lse = json.load(f)
#
#
#
# lse = LightStructureEnvironments.from_dict(dict_lse)
# pauling1=Pauling1_general_limit_rule(lse)
# print(pauling1.get_details())
# # # print(lse.structure)
# #
# # connection = Graph_Connection()
# # connection.newsetup(lse=lse)

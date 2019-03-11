from Graph_Test import Graph_Connection
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import json

# "mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-5986.json", "mp-306.json"
with open("tests/mp-7000.json", "r") as f:
    dict_lse = json.load(f)

lse = LightStructureEnvironments.from_dict(dict_lse)
# print(lse.structure)

connection = Graph_Connection()
connection.newsetup(lse=lse)

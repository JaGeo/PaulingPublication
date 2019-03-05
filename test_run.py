from PaulingRules import Pauling5
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import json

# "mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-5986.json", "mp-306.json"
with open("mp-19359.json", "r") as f:
    dict_lse = json.load(f)

lse = LightStructureEnvironments.from_dict(dict_lse)
# print(lse.structure)

pauling5 = Pauling5()
pauling5.newsetup(lse=lse, save_to_file=False)

print(pauling5.is_fulfilled())
#print(pauling5.get_details(options='CN'))
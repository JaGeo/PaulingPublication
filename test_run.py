from PaulingRules import Pauling3
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import json




#"mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-5986.json", "mp-306.json"
with open("mp-7000.json","r") as f:
    dict_lse=json.load(f)


lse=LightStructureEnvironments.from_dict(dict_lse)
#print(lse.structure)

pauling3=Pauling3(lse=lse,save_to_file=True,filename="mp-7000")

#schoenes ausgabeformat überlegen
print(pauling3.get_details())

print(pauling3.is_fulfilled())

#write a nice unit test - running through nearly all parts of the code



#brauche lse, um zu testen


from PaulingRules import Pauling4
from pymatgen.analysis.chemenv.coordination_environments.structure_environments import LightStructureEnvironments
import json

# "mp-1788.json", "mp-7000.json", "mp-19359.json", "mp-5986.json", "mp-306.json"
with open("mp-19359.json", "r") as f:
    dict_lse = json.load(f)

lse = LightStructureEnvironments.from_dict(dict_lse)
# print(lse.structure)

pauling4 = Pauling4()
pauling4.newsetup(lse=lse, save_to_file=False)

#print(pauling4.is_fulfilled())
#print(pauling4.get_details())
print(pauling4._postevaluation4thrule())
print(pauling4._postevaluation4thruleperpolyhedron_only_withoutproduct(6,6,1,3))
print(pauling4._postevaluation4thruleperpolyhedron_only_withoutproduct(6,6,3,1))
print(pauling4._postevaluation4thruleperpolyhedron_only_withoutproduct(6,6,1,1))
print(pauling4._postevaluation4thruleperpolyhedron_only_withoutproduct(6,6,3,3))
#teste from file and to file again
exit()
# schoenes ausgabeformat überlegen
#print(pauling3.get_details(maximumCN=4))

#print(pauling3.is_fulfilled())

# write a nice unit test - running through nearly all parts of the code


# brauche lse, um zu testen

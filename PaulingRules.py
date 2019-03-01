import json
from pymatgen.core.periodic_table import Element
from pymatgen.core.periodic_table import Specie
from pymatgen.core import PeriodicSite
import pymongo
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.magnetism.jahnteller import JahnTellerAnalyzer
from pymatgen.core.structure import Structure


# besserer aufbau der klassen?
# init: enthaelt lse und wenige spezifikationen
# get_details()
# is_fulfilled()

class Pauling1:
    # TODO: simplify class

    def __init__(self, lse, filenameradii="univalent_cat_radii.json", onlylowerlimit=False):
        """
        class to test Pauling's first rule
        :param lse: LightStructureEnvironment
        :param filenameradii: name of the file containing radius ratios
        :param onlylowerlimit: If False, ratio windows are considered
        """
        # important parameters
        self.cat_list = {}
        self.cat_valence_list = {}

        self.mat_pauling_fulfilled = 0
        self.env_not = 0
        self.no_env = 0
        self.no_cat = 0

        with open(filenameradii) as dd:
            dict_radii = json.load(dd)
        if self._checkonlyoxygenanion(lse):

            # TODO: could think of sorting lse.coordination_environments
            for isite, site_envs in enumerate(lse.coordination_environments):
                if site_envs != None:
                    if len(site_envs) > 0:
                        try:
                            try:
                                iratio = dict_radii[lse.structure.species[isite].symbol]['ratio_oxide']
                            except:
                                raise ValueError(
                                    lse.structure.species[isite].symbol + " not in the Pauling list")
                            if self.first_rule(iratio, site_envs[0], onlylowerlimit=onlylowerlimit):
                                cat_symbol = lse.structure.species[isite].symbol
                                val = lse.valences[isite]
                                # hier muss was rein
                                if not cat_symbol in self.cat_list:
                                    self.cat_list[cat_symbol] = [0, 0]
                                if not cat_symbol in self.cat_valence_list:
                                    self.cat_valence_list[cat_symbol] = {}
                                if not val in self.cat_valence_list[cat_symbol]:
                                    self.cat_valence_list[cat_symbol][val] = [0, 0]

                                self.cat_list[cat_symbol][0] = \
                                    self.cat_list[cat_symbol][0] + 1
                                self.cat_valence_list[cat_symbol][lse.valences[isite]][0] = \
                                    self.cat_valence_list[cat_symbol][lse.valences[isite]][0] + 1
                                self.mat_pauling_fulfilled += 1
                            else:
                                cat_symbol = lse.structure.species[isite].symbol
                                val = lse.valences[isite]
                                # hier muss was rein
                                if not cat_symbol in self.cat_list:
                                    self.cat_list[cat_symbol] = [0, 0]
                                if not cat_symbol in self.cat_valence_list:
                                    self.cat_valence_list[cat_symbol] = {}
                                if not val in self.cat_valence_list[cat_symbol]:
                                    self.cat_valence_list[cat_symbol][val] = [0, 0]

                                self.cat_list[lse.structure.species[isite].symbol][1] = \
                                    self.cat_list[lse.structure.species[isite].symbol][1] + 1
                                self.cat_valence_list[lse.structure.species[isite].symbol][lse.valences[isite]][1] = \
                                    self.cat_valence_list[lse.structure.species[isite].symbol][lse.valences[isite]][
                                        1] + 1

                                self.env_not += 1
                        except ValueError as err:
                            if err.args[0] == "env not in Pauling book":
                                self.no_env += 1
                            else:
                                self.no_cat += 1

    def get_details(self):
        """

        :return: returns a dict with many details on the Pauling rules
        """
        Outputdict = {}
        Outputdict["cat_dependency"] = self.cat_list
        Outputdict["cat_val_dependency"] = self.cat_valence_list
        Outputdict["Env_fulfilled"] = self.mat_pauling_fulfilled
        Outputdict["Env_notfulfilled"] = self.env_not
        Outputdict["Env_out_of_list"] = self.no_env
        Outputdict["Cat_out_of_list"] = self.no_cat

        return Outputdict

    def is_fulfilled(self):
        """
        tells you if the rule is fulfilled.
        :return: Boolean
        raises TypeError if not all environments and cations can be considered
        """
        if self.no_env != 0 or self.no_cat != 0:
            raise TypeError

        if self.env_not == 0 and self.mat_pauling_fulfilled > 0:
            return True
        else:
            return False

    # directly from Pauling's book
    def _predict_env_pauling_window(self, ratio):
        if ratio < 0.225:
            return ['does not exist']
        elif ratio >= 0.225 and ratio < 0.414:
            return ['T:4']  # tetrahedron
        elif ratio >= 0.414 and ratio < 0.592:
            return ['O:6']  # octahedron
        elif ratio >= 0.592 and ratio < 0.645:
            return ['FO:7']  # face-capped octahedron
        elif ratio >= 0.645 and ratio < 0.732:
            return ['SA:8']  # square antiprism
        elif ratio >= 0.732 and ratio < 1.0:
            return ['TT_1:9', 'C:8']  # Tricapped triangular prism (three square-face caps), #cube
        elif ratio >= 1.0:
            return ['C:12']  # Cuboctahedron

    def _predict_env_pauling_lowerlimit(self, ratio):
        ListToReturn = []
        if ratio >= 0.225:
            ListToReturn.append('T:4')  # tetrahedron
        if ratio >= 0.414:
            ListToReturn.append('O:6')  # octahedron
        if ratio >= 0.592:
            ListToReturn.append('FO:7')  # face-capped octahedron
        if ratio >= 0.645:
            ListToReturn.append('SA:8')  # square antiprism
        if ratio >= 0.732:
            ListToReturn.append('TT_1:9')  # Tricapped triangular prism (three square-face caps), #cube
            ListToReturn.append('C:8')
        if ratio >= 1.0:
            ListToReturn.append('C:12')  # Cuboctahedron
        return ListToReturn

    def _checkonlyoxygenanion(self, lse):
        """
        check if oxide
        :param lse: LightStructureEnvironment Object
        :return:
        """
        onlyoxygen = True
        sites = lse.structure.sites
        valences = lse.valences
        try:
            for isite, site in enumerate(sites):
                if sites[isite].species_string != 'O' and valences[isite] < 0:
                    onlyoxygen = False

        except:
            onlyoxygen = False
            print('Excluded material')
        return onlyoxygen

    def first_rule(self, iratio, site_env, onlylowerlimit):

        environments = ['T:4', 'O:6', 'FO:7',
                        'SA:8', 'TT_1:9', 'C:8', 'C:12']

        if site_env['ce_symbol'] in environments:
            # only values from Pauling book are considered
            if not onlylowerlimit:
                if str(site_env['ce_symbol']) in self._predict_env_pauling_window(iratio):
                    return True
                else:
                    return False

            else:
                if str(site_env['ce_symbol']) in self._predict_env_pauling_lowerlimit(iratio):
                    return True
                else:
                    return False

        else:
            raise ValueError("env not in Pauling book")

    # predict_env_pauling_lowerlimit_extended


# TODO: build in is_fulfilled, get_details

class Pauling2:
    def __init__(self, lse):
        """
        Class to test the electrostatic valence rule
        :param lse: LightStructureEnvironment
        """
        self.electrostatic_bond_strengths = {}
        for isite, site in enumerate(lse.structure):
            if lse.valences[isite] > 0:
                # print(site)
                # print(lse.coordination_environments[isite])
                nb_set = lse.neighbors_sets[isite][0]
                cn = float(len(nb_set))
                for nb_dict in nb_set.neighb_sites_and_indices:
                    # print(nb_dict)
                    nb_isite = nb_dict['index']
                    if nb_isite not in self.electrostatic_bond_strengths:
                        self.electrostatic_bond_strengths[nb_isite] = []
                    self.electrostatic_bond_strengths[nb_isite].append({'cation_isite': isite,
                                                                        'bond_strength': float(
                                                                            lse.valences[isite]) / cn})
        self.anions_bond_strengths = []
        self.satisfied = True
        self.lse = lse
        for isite, site in enumerate(lse.structure):
            if lse.valences[isite] < 0:
                if isite in self.electrostatic_bond_strengths:
                    bond_strengths = [nb['bond_strength']
                                      for nb in self.electrostatic_bond_strengths[isite]]
                    cations_isites = [nb['cation_isite']
                                      for nb in self.electrostatic_bond_strengths[isite]]
                else:
                    bond_strengths = []
                    cations_isites = []
                bond_strengths_sum = sum(bond_strengths)
                site_is_satisfied = bool(np.isclose(
                    bond_strengths_sum, -lse.valences[isite]))
                self.anions_bond_strengths.append({'anion_isite': isite,
                                                   'bond_strengths': bond_strengths,
                                                   'cations_isites': cations_isites,
                                                   'nominal_oxidation_state': lse.valences[isite],
                                                   'bond_strengths_sum': bond_strengths_sum,
                                                   'electrostatic_rule_satisfied': site_is_satisfied})
                if not site_is_satisfied:
                    self.satisfied = False

    def is_fulfilled(self):
        """
        Tells you if rule is fulfilled for the whole structure
        :return: Boolean
        """
        return self.satisfied

    def get_details(self):
        """
        Gives you an output dict with information on each anion
        :return: OutputDict with information on each anion
        """
        OutputDict = {}
        OutputDict["bvs_for_each_anion"] = self.get_anions_bvs()
        OutputDict["cations_around_anion"] = self.get_cations_around_anion()
        return OutputDict

    def get_anions_bvs(self):
        """
        get bond valences sums for each anion
        :return: list of bvs
        """
        bvs = []
        ianionsite = 0
        for isite, site in enumerate(self.lse.structure):
            if self.lse.valences[isite] < 0:
                bvs.append(
                    self.anions_bond_strengths[ianionsite]['bond_strengths_sum'])
                ianionsite = ianionsite + 1
        return bvs

    def get_cations_around_anion(self):
        """
        returns list of list with elements around anion
        :return:
        """
        ianionsite = 0
        elements = []
        for isite, site in enumerate(self.lse.structure):
            if self.lse.valences[isite] < 0:
                elements_around = [self.lse.structure[icat].species_string for icat in
                                   self.anions_bond_strengths[ianionsite]['cations_isites']]
                elements.append(elements_around)
                ianionsite = ianionsite + 1
        return elements


class PaulingConnection():
    """Class that can analyze connections of polyhedra"""

    def __init__(self, DISTANCE):
        self.DISTANCE = DISTANCE

    def _is_cationic_site(self, isite, valences):
        """

        :param isite: number of the site that is relevant
        :param valences: list of valences in order of sites
        :return:
        """
        if valences[isite] >= 0:
            return True
        else:
            return False

    def _get_oxygen_neighbors(self, lse, site, r, CN, ceindex):
        """
        :param lse: LightStructureEnvironment
        :param site: relevant site
        :param r: radius that should be considered for the analysis of the neighbors
        :param CN: coordination number of this cationic site
        :param ceindex: index of the coordination number at hand (for the most important one, usually 0)
        :return: list of oxygen neighbors
        """
        struct = lse.structure
        all_neighbors = struct.get_neighbors(site=site, r=r)
        sites_and_indices = lse.neighbors_sets[self._get_site_index(
            site, struct)][ceindex].neighb_sites_and_indices
        indices = [x['index'] for x in sites_and_indices]
        sorted_all = sorted(all_neighbors, key=lambda x: x[1])
        counter = 0
        oxygen_neighbors = []
        for i in sorted_all:
            if self._get_site_index(i[0], struct) in indices and counter < CN:
                oxygen_neighbors.append(i[0])
                counter = counter + 1
                indices.remove(self._get_site_index(i[0], struct))
            if counter == CN:
                break
        return oxygen_neighbors

    def _get_cation_neighbors(self, struct, site, r, valences=[]):
        """
        get cationic neighbors of a cation
        :param struct: Structure Object
        :param site: relevant cationic site
        :param r: radius
        :param valences: list of valences
        :return: list of neighboring sites
        """
        all_neighbors = struct.get_neighbors(
            site=site, r=r, include_index=False)
        cation_neighbors = []
        for i in all_neighbors:
            if self._is_cationic_site(self._get_site_index(i[0], struct), valences):
                cation_neighbors.append(i[0])
        return cation_neighbors

    def _get_site_index(self, insite, struct):
        sites = struct.sites
        for isite, site in enumerate(sites):
            if insite.is_periodic_image(site):
                return isite


class Pauling3and4(PaulingConnection):
    def __init__(self):
        pass
        # get distance by a parameter
        # self._newsetup(lse=lse, mat=filename, save=save_to_file, foldername=foldername)

    def newsetup(self, lse, filename, save_to_file=True, foldername='ThirdRuleAnalysisConnections', distance=8.0):
        """
            :param lse: LightStructureEnvironment
            :param save_to_file: Boolean
            :param filename: beginning of the file name, without ".json"
            :param foldername: name of the folder
            :param distance: float giving the distances of cations that is considered
        """
        super().__init__(DISTANCE=distance)
        if save_to_file:
            if not os.path.isdir(foldername):
                os.mkdir(foldername)
        DISTANCE = self.DISTANCE
        self.mat = filename
        lsegood = self._test_lse(lse)
        if lsegood:
            # get all relevant sites
            allsites = self._get_allsites(lse, DISTANCE)
            # search all pairs of sites that fulfill the DISTANCE criterion
            pairedsites = self._get_pairedsites(allsites, DISTANCE)
            self.PolyhedronDict = self._get_connections(pairedsites, lse, DISTANCE)
            if save_to_file:
                with open(os.join(foldername, self.mat + '.json'), 'w') as file:
                    json.dump(self.PolyhedronDict, file)

    def fromfile(self, filename, foldername='ThirdRuleAnalysisConnections'):
        # write it differently
        self.mat = filename
        with open(foldername + '/' + mat + '.json') as file:
            self.PolyhedronDict = json.load(file)

    def _get_number_connected_oxygens_and_CN(self, index, index2, sitetupel, lse):
        """

        :param index: site index of first cation site
        :param index2: ite index of second cation site
        :param sitetupel: tupel of sites that is investigated
        :param lse: LightStructureEnvironment
        :return: returns number of connected oxygens, coordination number of first cation and second cation
        """
        maxce = max(
            lse.coordination_environments[index], key=lambda x: x['ce_fraction'])
        numberce = lse.coordination_environments[index].index(maxce)
        CN = int(maxce['ce_symbol'].split(":")[1])
        oxygensites = self._get_oxygen_neighbors(
            lse, sitetupel[0], 8, CN, numberce)

        maxce2 = max(
            lse.coordination_environments[index2], key=lambda x: x['ce_fraction'])
        numberce2 = lse.coordination_environments[index].index(maxce)
        CN2 = int(maxce2['ce_symbol'].split(":")[1])

        oxygensites2 = self._get_oxygen_neighbors(
            lse, sitetupel[1], 8, CN2, numberce2)
        intersectOx = list(set(oxygensites).intersection(oxygensites2))
        numberconnect = len(intersectOx)
        return numberconnect, CN, CN2

    def _get_cationvalences(self, valences):
        """
        get valences of cation
        :param valences: valences of all sites
        :return: returns list of valences of cationic sites
        """
        cationvalences = []
        for i, val in enumerate(valences):
            if self._is_cationic_site(i, valences):
                cationvalences.append(val)
        return cationvalences

    def _get_cationCN(self, lse):
        """
        :param lse: LightStructureEnvironment
        :return: returns coordination numbers of cations as a list
        """
        sites = lse.structure.sites
        CNlist = []
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=lse.valences):
                maxce = max(
                    lse.coordination_environments[isite], key=lambda x: x['ce_fraction'])
                CN = int(maxce['ce_symbol'].split(":")[1])
                CNlist.append(CN)
        return CNlist

    def _get_connections(self, pairedsites, lse, DISTANCE):
        """

        :param pairedsites: paired sites
        :param lse: LightStructureEnvironments
        :param DISTANCE: distance for the connetions
        :return: Outputdict with many information
        """
        struct = lse.structure

        valences = []
        try:
            valences = lse.valences
            # print(valences)
        except:
            print('Valences not working')

        # get list of cation valences
        cationvalences = self._get_cationvalences(valences)
        # is the valence equal
        samevalences = (cationvalences.count(
            cationvalences[0]) == len(cationvalences))
        # get list of coordination numbers for each cation
        CNlist = self._get_cationCN(lse)
        sameCN = (CNlist.count(CNlist[0]) == len(CNlist))

        # do something else
        additionalinfo = []
        herepolyhedra = []
        for sitetupel in pairedsites:
            mydict = {}
            mydict['valences'] = []
            mydict['CN'] = []
            mydict['cations'] = []
            index = self._get_site_index(sitetupel[0], struct)
            index2 = self._get_site_index(sitetupel[1], struct)
            numberconnect, CN, CN2 = self._get_number_connected_oxygens_and_CN(
                index, index2, sitetupel, lse)
            try:
                mydict['valences'] = [valences[index], valences[index2]]
            except:
                pass
            mydict['CN'] = [CN, CN2]

            mydict['distance'] = np.linalg.norm(
                sitetupel[0].coords - sitetupel[1].coords)
            cat1 = sitetupel[0].species_string
            cat2 = sitetupel[1].species_string
            mydict['cations'] = [cat1, cat2]
            additionalinfo.append(mydict)
            herepolyhedra.append(numberconnect)

        numberofPolyhedra = len(herepolyhedra)
        nofPoly_notconnected = len([x for x in herepolyhedra if x == 0])
        nofPoly_corner = len([x for x in herepolyhedra if x == 1])
        noofPoly_edge = len([x for x in herepolyhedra if x == 2])
        noofPoly_face = len([x for x in herepolyhedra if x >= 3])

        outputdict = {}
        outputdict['MaxDistance'] = DISTANCE
        outputdict['Polyhedra'] = numberofPolyhedra
        outputdict['Not'] = nofPoly_notconnected
        outputdict['Corner'] = nofPoly_corner
        outputdict['Edge'] = noofPoly_edge
        outputdict['Face'] = noofPoly_face
        outputdict['PolyConnect'] = herepolyhedra
        outputdict['Additional'] = additionalinfo
        outputdict['cationvalences'] = cationvalences
        outputdict['samevalences'] = samevalences
        outputdict['CNlist'] = CNlist
        outputdict['sameCN'] = sameCN

        return outputdict

    def _test_lse(self, lse):
        """
        tests whether a lse is present for each cation
        :param lse: LightStructureEnvironment
        :return: Boolean
        """
        lsegood = True
        for isite, site in enumerate(lse.structure.sites):
            if self._is_cationic_site(isite, valences=lse.valences):
                try:
                    lse.coordination_environments[isite]
                except:
                    print('Excluded:' + mat)
                    lsegood = False
                    break
        return lsegood

    def _get_pairedsites(self, allsites, DISTANCE):
        """

        :param allsites: all sites that will be considered
        :param DISTANCE: distance that is considered
        :return: paired sites with a certain distance
        """
        startendpoints = []
        pairedsites = []
        for i in range(0, len(allsites)):
            for j in range(i, len(allsites)):
                #
                fraccoorda = allsites[i].frac_coords
                fraccoordb = allsites[j].frac_coords

                if not (np.linalg.norm(allsites[i].coords - allsites[j].coords) > DISTANCE or (
                        fraccoorda[0] == fraccoordb[0] and fraccoorda[1] == fraccoordb[1] and fraccoorda[2] ==
                        fraccoordb[2])):

                    # make sure that at least one coordinate lies in the first cell [0:1)
                    while fraccoorda[0] < 0.0 or fraccoorda[1] < 0.0 or fraccoorda[2] < 0.0 or fraccoordb[
                        0] < 0.0 or fraccoordb[1] < 0.0 or fraccoordb[2] < 0.0:
                        # print('test')
                        for ij in range(0, 3):
                            if fraccoorda[ij] < 0.0 or fraccoordb[ij] < 0.0:
                                fraccoorda[ij] = fraccoorda[ij] + 1.0
                                fraccoordb[ij] = fraccoordb[ij] + 1.0
                    while not ((fraccoorda[0] >= 0.0 and fraccoordb[0] >= 0.0)
                               and (fraccoorda[0] < 1.0 or fraccoordb[0] < 1.0)
                               and (fraccoorda[1] >= 0.0 and fraccoordb[1] >= 0.0)
                               and (fraccoorda[1] < 1.0 or fraccoordb[1] < 1.0)
                               and (fraccoorda[2] >= 0.0 and fraccoordb[2] >= 0.0)
                               and (fraccoorda[2] < 1.0 or fraccoordb[2] < 1.0)):
                        for number in range(0, 3):
                            if ((fraccoorda[number] >= 1.0) and fraccoordb[number] >= 1.0):
                                # print('test2')
                                fraccoorda[number] = fraccoorda[number] - 1.0
                                fraccoordb[number] = fraccoordb[number] - 1.0

                    # Sort to decide whether we already considered this pair
                    sort = sorted([fraccoorda, fraccoordb],
                                  key=lambda x: (x[0], x[1], x[2]))
                    contained = False
                    for istartendpoints in startendpoints:
                        if np.allclose(sort, istartendpoints):
                            contained = True
                            break
                    if not contained:
                        startendpoints.append(sort)
                        pairedsites.append([allsites[i], allsites[j]])
        return pairedsites

    def _get_allsites(self, lse, DISTANCE):
        """

        :param lse: LightStructureEnvironment
        :param DISTANCE: distance between cations
        :return: returns relevant sites as a list
        """

        sites = lse.structure.sites
        valences = lse.valences
        struct = lse.structure
        allsites = []
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=valences):
                # TODO: reprogram this part!
                # get infos from site
                atoms_n_occu = site.species_string
                lattice = site.lattice
                coords = site.frac_coords

                # get number of cells that have to be considered to find neighbors with a CERTAIN distance
                ceila = int(np.ceil(DISTANCE / struct.lattice.a))
                ceilb = int(np.ceil(DISTANCE / struct.lattice.b))
                ceilc = int(np.ceil(DISTANCE / struct.lattice.c))

                # +1 fore safety
                # TODO: take one with
                # could be done in a safer way
                # TODO: do it safer!!!
                for i in range(0, (ceila + 1)):
                    for j in range(0, (ceilb + 1)):
                        for k in range(0, ceilc + 1):
                            addcoord = [float(i), float(j), float(k)]
                            newcoords = coords + addcoord
                            allsites.append(PeriodicSite(
                                atoms_n_occu, newcoords, lattice, to_unit_cell=False))

        return allsites

    def get_connections(self):
        """
        Gives connections
        :return: Outputdict with number of connections
        """
        OutputDict={}
        OutputDict['no']=self.PolyhedronDict['Not']
        OutputDict['corner']=self.PolyhedronDict['Corner']
        OutputDict['edge']=self.PolyhedronDict['Edge']
        OutputDict['face']=self.PolyhedronDict['Face']
        return OutputDict


class Pauling3(Pauling3and4):
    def is_fulfilled(self):
        # TODO: Kritrien, die erfuellt sein muessen
        # only connections via corner and edges, no face connections
        pass

    def get_details(self):
        # get details on the evaluation
        pass

    def postevaluation3rdrule(self, maxCN=None):
        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        allPolyhedra = len(herepolyhedra)
        Polyhedra_notconnected = len([x for x in herepolyhedra if x == 0])
        Polyhedra_corner = len([x for x in herepolyhedra if x == 1])
        Polyhedra_edge = len([x for x in herepolyhedra if x == 2])
        Polyhedra_face = len([x for x in herepolyhedra if x >= 3])

        nocoulomb = []
        cornercoulomb = []
        edgecoulomb = []
        facecoulomb = []
        nodistances = []
        cornerdistances = []
        edgedistances = []
        facedistances = []

        highvalence_notconnected = 0
        highvalence_corner = 0
        highvalence_edge = 0
        highvalence_face = 0
        lowvalence_notconnected = 0
        lowvalence_corner = 0
        lowvalence_edge = 0
        lowvalence_face = 0

        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if ((iinfo['valences'][0] * iinfo['valences'][1] > 4) and (iinfo['CN'][0] * iinfo['CN'][1] < 36)):

                if herepolyhedra[numberpolyhedra] == 0:
                    highvalence_notconnected = highvalence_notconnected + 1
                    nocoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    nodistances.append(iinfo['distance'])
                elif herepolyhedra[numberpolyhedra] == 1:
                    highvalence_corner = highvalence_corner + 1
                    cornercoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    cornerdistances.append(iinfo['distance'])
                elif herepolyhedra[numberpolyhedra] == 2:
                    highvalence_edge = highvalence_edge + 1
                    edgecoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    edgedistances.append(iinfo['distance'])
                elif herepolyhedra[numberpolyhedra] > 2:
                    highvalence_face = highvalence_face + 1
                    facecoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    facedistances.append(iinfo['distance'])
            else:
                if herepolyhedra[numberpolyhedra] == 0:
                    lowvalence_notconnected = lowvalence_notconnected + 1
                    nocoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    nodistances.append(iinfo['distance'])
                elif herepolyhedra[numberpolyhedra] == 1:
                    lowvalence_corner = lowvalence_corner + 1
                    cornercoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    cornerdistances.append(iinfo['distance'])
                elif herepolyhedra[numberpolyhedra] == 2:
                    lowvalence_edge = lowvalence_edge + 1
                    edgecoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    edgedistances.append(iinfo['distance'])
                elif herepolyhedra[numberpolyhedra] > 2:
                    lowvalence_face = lowvalence_face + 1
                    facecoulomb.append(
                        iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                    facedistances.append(iinfo['distance'])
            numberpolyhedra = numberpolyhedra + 1

        if maxCN is not None:
            max_CN_all = 0
            maxCN_notconnected = 0
            maxCN_corner = 0
            maxCN_edge = 0
            maxCN_face = 0
            numberpolyhedra = 0
            for iinfo in additionalinfo:
                if (iinfo['CN'][0] <= 8 and iinfo['CN'][1] <= 8):

                    max_CN_all = max_CN_all + 1

                    if herepolyhedra[numberpolyhedra] == 0:
                        maxCN_notconnected = maxCN_notconnected + 1

                    elif herepolyhedra[numberpolyhedra] == 1:
                        maxCN_corner = maxCN_corner + 1

                    elif herepolyhedra[numberpolyhedra] == 2:
                        maxCN_edge = maxCN_edge + 1

                    elif herepolyhedra[numberpolyhedra] > 2:
                        maxCN_face = maxCN_face + 1

                numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['Polyhedra Pairs'] = allPolyhedra
        Thirdruledict['not connected'] = Polyhedra_notconnected
        Thirdruledict['corner sharing'] = Polyhedra_corner
        Thirdruledict['edge sharing'] = Polyhedra_edge
        Thirdruledict['face sharing'] = Polyhedra_face
        Thirdruledict['high valence not connected'] = highvalence_notconnected
        Thirdruledict['high valence corner sharing'] = highvalence_corner
        Thirdruledict['high valence edge sharing'] = highvalence_edge
        Thirdruledict['high valence face sharing'] = highvalence_face
        Thirdruledict['low valence not connected'] = lowvalence_notconnected
        Thirdruledict['low valence corner sharing'] = lowvalence_corner
        Thirdruledict['low valence edge sharing'] = lowvalence_edge
        Thirdruledict['low valence face sharing'] = lowvalence_face
        Thirdruledict['q1q2/r12 not connected'] = nocoulomb
        Thirdruledict['q1q2/r12 corner sharing'] = cornercoulomb
        Thirdruledict['q1q2/r12 edge sharing'] = edgecoulomb
        Thirdruledict['q1q2/r12 face sharing'] = facecoulomb
        Thirdruledict['distances not connected'] = nodistances
        Thirdruledict['distances corner sharing'] = cornerdistances
        Thirdruledict['distances edge sharing'] = edgedistances
        Thirdruledict['distances face sharing'] = facedistances

        if maxCN is not None:
            Thirdruledict['Polyhedra Pairs maxCN'] = max_CN_all
            Thirdruledict['not connected maxCN'] = maxCN_notconnected
            Thirdruledict['corner sharing maxCN'] = maxCN_corner
            Thirdruledict['edge sharing maxCN'] = maxCN_edge
            Thirdruledict['face sharing maxCN'] = maxCN_face

        return Thirdruledict

    def postevaluation3rdrule_valenceandCNdependend(self, LowerLimitValence=4, UpperLimitCN=36):
        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        highvalence_notconnected = 0
        highvalence_corner = 0
        highvalence_edge = 0
        highvalence_face = 0
        lowvalence_notconnected = 0
        lowvalence_corner = 0
        lowvalence_edge = 0
        lowvalence_face = 0

        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if ((iinfo['valences'][0] * iinfo['valences'][1] > LowerLimitValence) and (
                    iinfo['CN'][0] * iinfo['CN'][1] < UpperLimitCN)):

                if herepolyhedra[numberpolyhedra] == 0:
                    highvalence_notconnected = highvalence_notconnected + 1

                elif herepolyhedra[numberpolyhedra] == 1:
                    highvalence_corner = highvalence_corner + 1

                elif herepolyhedra[numberpolyhedra] == 2:
                    highvalence_edge = highvalence_edge + 1

                elif herepolyhedra[numberpolyhedra] > 2:
                    highvalence_face = highvalence_face + 1

            else:
                if herepolyhedra[numberpolyhedra] == 0:
                    lowvalence_notconnected = lowvalence_notconnected + 1

                elif herepolyhedra[numberpolyhedra] == 1:
                    lowvalence_corner = lowvalence_corner + 1

                elif herepolyhedra[numberpolyhedra] == 2:
                    lowvalence_edge = lowvalence_edge + 1

                elif herepolyhedra[numberpolyhedra] > 2:
                    lowvalence_face = lowvalence_face + 1

            numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['high valence not connected'] = highvalence_notconnected
        Thirdruledict['high valence corner sharing'] = highvalence_corner
        Thirdruledict['high valence edge sharing'] = highvalence_edge
        Thirdruledict['high valence face sharing'] = highvalence_face
        Thirdruledict['low valence not connected'] = lowvalence_notconnected
        Thirdruledict['low valence corner sharing'] = lowvalence_corner
        Thirdruledict['low valence edge sharing'] = lowvalence_edge
        Thirdruledict['low valence face sharing'] = lowvalence_face

        return Thirdruledict

    def postevaluation3rdrule_valenceandCNdependend_nogroups(self, ValenceProduct, CNProduct, maxCN=None):
        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        notconnected = 0
        corner = 0
        edge = 0
        face = 0

        # Dict that includes cations sharing no connections, corners, edges, and faces
        # This can be used to evaluate the results per polyhedron pair
        Dict_ElementDependency = {}

        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if maxCN is None:
                if ((iinfo['valences'][0] * iinfo['valences'][1] == ValenceProduct) and (
                        iinfo['CN'][0] * iinfo['CN'][1] == CNProduct)):
                    # if iinfo['cations'][0] not in Dict_ElementDependency:
                    #     Dict_ElementDependency[iinfo['cations'][0]]={'no':0,'corner':0,'edge':0,'face':0}
                    #
                    # if iinfo['cations'][1] not in Dict_ElementDependency:
                    #     Dict_ElementDependency[iinfo['cations'][1]]={'no':0,'corner':0,'edge':0,'face':0}

                    if herepolyhedra[numberpolyhedra] == 0:
                        notconnected = notconnected + 1
                        # Dict_ElementDependency[iinfo['cations'][0]]['no']+=1
                        # Dict_ElementDependency[iinfo['cations'][1]]['no']+=1

                    elif herepolyhedra[numberpolyhedra] == 1:
                        corner = corner + 1
                        # Dict_ElementDependency[iinfo['cations'][0]]['corner'] += 1
                        # Dict_ElementDependency[iinfo['cations'][1]]['corner'] += 1


                    elif herepolyhedra[numberpolyhedra] == 2:
                        edge = edge + 1
                        # Dict_ElementDependency[iinfo['cations'][0]]['edge'] += 1
                        # Dict_ElementDependency[iinfo['cations'][1]]['edge'] += 1


                    elif herepolyhedra[numberpolyhedra] > 2:
                        face = face + 1
                        # Dict_ElementDependency[iinfo['cations'][0]]['face'] += 1
                        # Dict_ElementDependency[iinfo['cations'][1]]['face'] += 1

            else:
                if iinfo["CN"][0] <= maxCN and iinfo["CN"][1] <= maxCN:
                    if ((iinfo['valences'][0] * iinfo['valences'][1] == ValenceProduct) and (
                            iinfo['CN'][0] * iinfo['CN'][1] == CNProduct)):
                        # if iinfo['cations'][0] not in Dict_ElementDependency:
                        #     Dict_ElementDependency[iinfo['cations'][0]] = {'no':0,'corner':0,'edge':0,'face':0}
                        # if iinfo['cations'][1] not in Dict_ElementDependency:
                        #     Dict_ElementDependency[iinfo['cations'][1]] = {'no':0,'corner':0,'edge':0,'face':0}

                        if herepolyhedra[numberpolyhedra] == 0:
                            notconnected = notconnected + 1
                            # Dict_ElementDependency[iinfo['cations'][0]]['no'] += 1
                            # Dict_ElementDependency[iinfo['cations'][1]]['no'] += 1

                        elif herepolyhedra[numberpolyhedra] == 1:
                            corner = corner + 1
                            # Dict_ElementDependency[iinfo['cations'][0]]['corner'] += 1
                            # Dict_ElementDependency[iinfo['cations'][1]]['corner'] += 1


                        elif herepolyhedra[numberpolyhedra] == 2:
                            edge = edge + 1
                            # Dict_ElementDependency[iinfo['cations'][0]]['edge'] += 1
                            # Dict_ElementDependency[iinfo['cations'][1]]['edge'] += 1


                        elif herepolyhedra[numberpolyhedra] > 2:
                            face = face + 1
                            # Dict_ElementDependency[iinfo['cations'][0]]['face'] += 1
                            # Dict_ElementDependency[iinfo['cations'][1]]['face'] += 1

            numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['not connected'] = notconnected
        Thirdruledict['corner sharing'] = corner
        Thirdruledict['edge sharing'] = edge
        Thirdruledict['face sharing'] = face
        # Thirdruledict['ElementDependency']=Dict_ElementDependency

        return Thirdruledict

    def postevaluation3rdrule_evaluate_elementdependency(self, maxCN=None):
        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        # notconnected = 0
        # corner = 0
        # edge = 0
        # face = 0

        # Dict that includes cations sharing no connections, corners, edges, and faces
        # This can be used to evaluate the results per polyhedron pair
        Dict_ElementDependency = {}

        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if maxCN is None:
                if iinfo['cations'][0] not in Dict_ElementDependency:
                    Dict_ElementDependency[iinfo['cations'][0]] = {'no': 0, 'corner': 0, 'edge': 0, 'face': 0}

                if iinfo['cations'][1] not in Dict_ElementDependency:
                    Dict_ElementDependency[iinfo['cations'][1]] = {'no': 0, 'corner': 0, 'edge': 0, 'face': 0}

                if herepolyhedra[numberpolyhedra] == 0:
                    Dict_ElementDependency[iinfo['cations'][0]]['no'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]]['no'] += 1

                elif herepolyhedra[numberpolyhedra] == 1:
                    Dict_ElementDependency[iinfo['cations'][0]]['corner'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]]['corner'] += 1


                elif herepolyhedra[numberpolyhedra] == 2:
                    Dict_ElementDependency[iinfo['cations'][0]]['edge'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]]['edge'] += 1


                elif herepolyhedra[numberpolyhedra] > 2:
                    Dict_ElementDependency[iinfo['cations'][0]]['face'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]]['face'] += 1

            else:
                if iinfo["CN"][0] <= maxCN and iinfo["CN"][1] <= maxCN:
                    if iinfo['cations'][0] not in Dict_ElementDependency:
                        Dict_ElementDependency[iinfo['cations'][0]] = {'no': 0, 'corner': 0, 'edge': 0, 'face': 0}
                    if iinfo['cations'][1] not in Dict_ElementDependency:
                        Dict_ElementDependency[iinfo['cations'][1]] = {'no': 0, 'corner': 0, 'edge': 0, 'face': 0}

                    if herepolyhedra[numberpolyhedra] == 0:
                        Dict_ElementDependency[iinfo['cations'][0]]['no'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]]['no'] += 1

                    elif herepolyhedra[numberpolyhedra] == 1:
                        Dict_ElementDependency[iinfo['cations'][0]]['corner'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]]['corner'] += 1


                    elif herepolyhedra[numberpolyhedra] == 2:
                        Dict_ElementDependency[iinfo['cations'][0]]['edge'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]]['edge'] += 1


                    elif herepolyhedra[numberpolyhedra] > 2:
                        Dict_ElementDependency[iinfo['cations'][0]]['face'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]]['face'] += 1

            numberpolyhedra = numberpolyhedra + 1

        return Dict_ElementDependency


class Pauling4(Pauling3and4):
    def postevaluation4thrule(self):
        inputdict = self.PolyhedronDict

        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']
        cationvalences = inputdict['cationvalences']
        samevalences = inputdict['samevalences']
        CNlist = inputdict['CNlist']
        sameCN = inputdict['sameCN']

        PaulingRuleCanBeTested = False

        #
        if (not samevalences) or (not sameCN):

            print(cationvalences)
            print(CNlist)
            maxCN = max(CNlist)
            minCN = min(CNlist)
            meanCN = (maxCN + minCN) / 2.0
            print(meanCN)
            maxVal = max(cationvalences)
            minVal = min(cationvalences)
            meanVal = (maxVal + minVal) / 2.0
            print(meanVal)
            smallnohere = 0
            smallcornerhere = 0
            smalledgehere = 0
            smallfacehere = 0
            numberpolyhedra = 0

            for iinfo in additionalinfo:
                if (iinfo['valences'][0] == maxVal and iinfo['CN'][0] == minCN) and (
                        iinfo['valences'][1] == maxVal and iinfo['CN'][1] == minCN):

                    if herepolyhedra[numberpolyhedra] == 0:
                        smallnohere = smallnohere + 1
                        print('no')
                    elif herepolyhedra[numberpolyhedra] == 1:
                        smallcornerhere = smallcornerhere + 1
                        print('corner')
                    elif herepolyhedra[numberpolyhedra] == 2:
                        smalledgehere = smalledgehere + 1
                        print('edge')
                    elif herepolyhedra[numberpolyhedra] > 2:
                        smallfacehere = smallfacehere + 1
                        print('face')

                numberpolyhedra = numberpolyhedra + 1

            # save number of structures fulfilling the rule( that means only those with high coordition number and small valences built connections)
            if smallcornerhere == 0 and smalledgehere == 0 and smallfacehere == 0:
                Pauling4thruleokay = True
                print('Pauling 4 okay')
            else:
                Pauling4thruleokay = False
                print('Pauling 4 not okay')
            PaulingRuleCanBeTested = True

        FourthRuleDict = {}
        if PaulingRuleCanBeTested:
            FourthRuleDict['Can be tested'] = True
            if Pauling4thruleokay:
                FourthRuleDict['Positive Test'] = True
            else:
                FourthRuleDict['Positive Test'] = False

        else:
            FourthRuleDict['Can be tested'] = False

        return FourthRuleDict

    def postevaluation4thruleperpolyhedron(self, CNlimit=36, valencelimit=4):
        # TODO: include CNmax

        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        allPolyhedra = len(herepolyhedra)
        Polyhedra_notconnected = len([x for x in herepolyhedra if x == 0])
        Polyhedra_corner = len([x for x in herepolyhedra if x == 1])
        Polyhedra_edge = len([x for x in herepolyhedra if x == 2])
        Polyhedra_face = len([x for x in herepolyhedra if x >= 3])

        samevalences = inputdict['samevalences']
        sameCN = inputdict['sameCN']

        #

        nocoulomb = []
        cornercoulomb = []
        edgecoulomb = []
        facecoulomb = []
        nodistances = []
        cornerdistances = []
        edgedistances = []
        facedistances = []

        highvalence_notconnected = 0
        highvalence_corner = 0
        highvalence_edge = 0
        highvalence_face = 0
        lowvalence_notconnected = 0
        lowvalence_corner = 0
        lowvalence_edge = 0
        lowvalence_face = 0

        numberpolyhedra = 0
        if (not samevalences) or (not sameCN):
            # correctly evaluated?

            for iinfo in additionalinfo:
                if ((iinfo['valences'][0] * iinfo['valences'][1] > valencelimit) and (
                        iinfo['CN'][0] * iinfo['CN'][1] < CNlimit)):

                    if herepolyhedra[numberpolyhedra] == 0:
                        highvalence_notconnected = highvalence_notconnected + 1
                        nocoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        nodistances.append(iinfo['distance'])
                    elif herepolyhedra[numberpolyhedra] == 1:
                        highvalence_corner = highvalence_corner + 1
                        cornercoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        cornerdistances.append(iinfo['distance'])
                    elif herepolyhedra[numberpolyhedra] == 2:
                        highvalence_edge = highvalence_edge + 1
                        edgecoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        edgedistances.append(iinfo['distance'])
                    elif herepolyhedra[numberpolyhedra] > 2:
                        highvalence_face = highvalence_face + 1
                        facecoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        facedistances.append(iinfo['distance'])
                else:
                    if herepolyhedra[numberpolyhedra] == 0:
                        lowvalence_notconnected = lowvalence_notconnected + 1
                        nocoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        nodistances.append(iinfo['distance'])
                    elif herepolyhedra[numberpolyhedra] == 1:
                        lowvalence_corner = lowvalence_corner + 1
                        cornercoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        cornerdistances.append(iinfo['distance'])
                    elif herepolyhedra[numberpolyhedra] == 2:
                        lowvalence_edge = lowvalence_edge + 1
                        edgecoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        edgedistances.append(iinfo['distance'])
                    elif herepolyhedra[numberpolyhedra] > 2:
                        lowvalence_face = lowvalence_face + 1
                        facecoulomb.append(
                            iinfo['valences'][0] * iinfo['valences'][1] / iinfo['distance'])
                        facedistances.append(iinfo['distance'])
                numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['Polyhedra Pairs'] = allPolyhedra
        Thirdruledict['not connected'] = Polyhedra_notconnected
        Thirdruledict['corner sharing'] = Polyhedra_corner
        Thirdruledict['edge sharing'] = Polyhedra_edge
        Thirdruledict['face sharing'] = Polyhedra_face
        Thirdruledict['high valence not connected'] = highvalence_notconnected
        Thirdruledict['high valence corner sharing'] = highvalence_corner
        Thirdruledict['high valence edge sharing'] = highvalence_edge
        Thirdruledict['high valence face sharing'] = highvalence_face
        Thirdruledict['low valence not connected'] = lowvalence_notconnected
        Thirdruledict['low valence corner sharing'] = lowvalence_corner
        Thirdruledict['low valence edge sharing'] = lowvalence_edge
        Thirdruledict['low valence face sharing'] = lowvalence_face
        Thirdruledict['q1q2/r12 not connected'] = nocoulomb
        Thirdruledict['q1q2/r12 corner sharing'] = cornercoulomb
        Thirdruledict['q1q2/r12 edge sharing'] = edgecoulomb
        Thirdruledict['q1q2/r12 face sharing'] = facecoulomb
        Thirdruledict['distances not connected'] = nodistances
        Thirdruledict['distances corner sharing'] = cornerdistances
        Thirdruledict['distances edge sharing'] = edgedistances
        Thirdruledict['distances face sharing'] = facedistances

        return Thirdruledict

    def postevaluation4thruleperpolyhedron_only(self, CNProduct, ValenceProduct):
        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        allPolyhedra = len(herepolyhedra)
        Polyhedra_notconnected = len([x for x in herepolyhedra if x == 0])
        Polyhedra_corner = len([x for x in herepolyhedra if x == 1])
        Polyhedra_edge = len([x for x in herepolyhedra if x == 2])
        Polyhedra_face = len([x for x in herepolyhedra if x >= 3])

        samevalences = inputdict['samevalences']
        sameCN = inputdict['sameCN']

        notconnected = 0
        corner = 0
        edge = 0
        face = 0

        numberpolyhedra = 0
        if (not samevalences) or (not sameCN):
            # is this really okay like this?
            # or should I evaluate not samevalences and not same CN seperately?

            for iinfo in additionalinfo:
                if ((iinfo['valences'][0] * iinfo['valences'][1] == ValenceProduct) and (
                        iinfo['CN'][0] * iinfo['CN'][1] == CNProduct)):

                    if herepolyhedra[numberpolyhedra] == 0:
                        notconnected = notconnected + 1
                    elif herepolyhedra[numberpolyhedra] == 1:
                        corner = corner + 1
                    elif herepolyhedra[numberpolyhedra] == 2:
                        edge = edge + 1
                    elif herepolyhedra[numberpolyhedra] > 2:
                        face = face + 1

                numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['Polyhedra Pairs'] = allPolyhedra
        Thirdruledict['not connected'] = Polyhedra_notconnected
        Thirdruledict['corner sharing'] = Polyhedra_corner
        Thirdruledict['edge sharing'] = Polyhedra_edge
        Thirdruledict['face sharing'] = Polyhedra_face
        Thirdruledict['not connected'] = notconnected
        Thirdruledict['corner sharing'] = corner
        Thirdruledict['edge sharing'] = edge
        Thirdruledict['face sharing'] = face

        return Thirdruledict

    # Add additional evaluation of 4th rule that does not depend on valence or CN product but on 2 valences, 2 CN
    def postevaluation4thruleperpolyhedron_only_withoutproduct(self, CN1, CN2, val1, val2):
        inputdict = self.PolyhedronDict
        if self.DISTANCE != inputdict['MaxDistance']:
            raise ValueError

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        allPolyhedra = len(herepolyhedra)
        Polyhedra_notconnected = len([x for x in herepolyhedra if x == 0])
        Polyhedra_corner = len([x for x in herepolyhedra if x == 1])
        Polyhedra_edge = len([x for x in herepolyhedra if x == 2])
        Polyhedra_face = len([x for x in herepolyhedra if x >= 3])

        samevalences = inputdict['samevalences']
        sameCN = inputdict['sameCN']

        notconnected = 0
        corner = 0
        edge = 0
        face = 0

        numberpolyhedra = 0
        if (not samevalences) or (not sameCN):
            # is this really okay like this?
            # or should I evaluate not samevalences and not same CN seperately?

            for iinfo in additionalinfo:
                # symmetry of valence
                if ((iinfo['valences'][0] == val1) and (iinfo['valences'][1] == val2)) or (
                        (iinfo['valences'][1] == val1) and (iinfo['valences'][0] == val2)):
                    if ((iinfo['CN'][0] == CN1 and iinfo['CN'][1] == CN2) or (
                            iinfo['CN'][0] == CN2 and iinfo['CN'][1] == CN1)):
                        if herepolyhedra[numberpolyhedra] == 0:
                            notconnected = notconnected + 1
                        elif herepolyhedra[numberpolyhedra] == 1:
                            corner = corner + 1
                        elif herepolyhedra[numberpolyhedra] == 2:
                            edge = edge + 1
                        elif herepolyhedra[numberpolyhedra] > 2:
                            face = face + 1

                        numberpolyhedra = numberpolyhedra + 1

        Thirdruledict = {}
        Thirdruledict['Polyhedra Pairs'] = allPolyhedra
        Thirdruledict['not connected'] = Polyhedra_notconnected
        Thirdruledict['corner sharing'] = Polyhedra_corner
        Thirdruledict['edge sharing'] = Polyhedra_edge
        Thirdruledict['face sharing'] = Polyhedra_face
        Thirdruledict['not connected'] = notconnected
        Thirdruledict['corner sharing'] = corner
        Thirdruledict['edge sharing'] = edge
        Thirdruledict['face sharing'] = face

        return Thirdruledict


class Pauling5(PaulingConnection):

    def __init(self, lse, filename, foldername):
        super().__init__(DISTANCE=distance)
        self.newsetup(lse, filename, foldername)

    def newsetup(self, lse, mat, foldername):
        # collects the coordination numbers, coordination environments and the number of connections via corners, edges and faces of each of the polyhedra for each of the cations

        struct = lse.structure
        sites = struct.sites

        catid = []
        catenv = []
        connection_corners = []
        connection_edges = []
        connection_faces = []
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=lse.valences):
                corner = 0
                edge = 0
                face = 0
                catid.append([site.species_string, lse.valences[isite]])
                maxce = max(
                    lse.coordination_environments[isite], key=lambda x: x['ce_fraction'])
                numberce = lse.coordination_environments[isite].index(maxce)
                CN = int(maxce['ce_symbol'].split(":")[1])
                oxygensites = self._get_oxygen_neighbors(
                    lse, site, 8, CN, numberce)
                catenv.append(maxce['ce_symbol'])
                cation_neighbors = self._get_cation_neighbors(
                    struct, site, r=self.DISTANCE, valences=lse.valences)
                for neigh in cation_neighbors:
                    indexneigh = self._get_site_index(neigh, struct)
                    maxce2 = max(
                        lse.coordination_environments[indexneigh], key=lambda x: x['ce_fraction'])
                    numberce2 = lse.coordination_environments[indexneigh].index(
                        maxce2)
                    CN2 = int(maxce2['ce_symbol'].split(":")[1])
                    oxygensites2 = self._get_oxygen_neighbors(lse, neigh, 8, CN2, numberce2)
                    intersectOx = list(
                        set(oxygensites).intersection(oxygensites2))
                    numberconnect = len(intersectOx)
                    if numberconnect == 1:
                        corner = corner + 1
                    if numberconnect == 2:
                        edge = edge + 1
                    if numberconnect >= 3:
                        face = face + 1
                connection_corners.append(corner)
                connection_edges.append(edge)
                connection_faces.append(face)

        # print(connections)
        # get unique catsite
        # erst mal schauen, ob es uberhaupt mehr als 1 Kation mit entsprechenden gleicher Elementsorte und CN gibt
        # check if there is more than one cation of the same element
        uniquecat = []
        for cat in catid:
            if cat not in uniquecat:
                uniquecat.append(cat)
        print(uniquecat)

        outputdict = {}
        outputdict['connection_corners'] = connection_corners
        outputdict['connection_edges'] = connection_edges
        outputdict['connection_faces'] = connection_faces
        outputdict['catenv'] = catenv
        outputdict['catid'] = catid
        outputdict['uniquecat'] = uniquecat

        with open(foldername + '/' + mat + '.json', 'w') as file:
            json.dump(outputdict, file)

        self.FifthRuleDict = outputdict

    def fromfile(self, mat, foldername):
        self.mat = mat
        print(mat)
        with open(foldername + '/' + mat + '.json') as file:
            self.FifthRuleDict = json.load(file)

    def postevaluation5thrule(self, mat, excluded, maxCN=None):
        # ToDO: build in max-CN
        connection_corners = self.FifthRuleDict['connection_corners']
        connection_edges = self.FifthRuleDict['connection_edges']
        connection_faces = self.FifthRuleDict['connection_faces']
        catenv = self.FifthRuleDict['catenv']
        catid = self.FifthRuleDict['catid']
        uniquecat = self.FifthRuleDict['uniquecat']

        hassamechemenvandsameconnectionnumber = True
        hassameenv = True
        hassameCN = True

        # print(catid)
        # print(uniquecat)

        skip = True

        # Checks if one of the unique cations is not in excluded and then proceeds
        # That's how it should be
        # Maybe make excluded also valence dependent.

        for unicat in uniquecat:
            if str(unicat) not in excluded:
                print(str(unicat) + ' is not in excluded')
                skip = False
            else:
                print(unicat)
        if skip:
            print(mat + " is excluded")

        if not skip:
            # tells you if all cations that are the same have the same CN, CE, CE with no connections
            if len(uniquecat) != len(catid):

                for icat, cat in enumerate(catid):
                    for icat2, cat2 in enumerate(catid):
                        if icat2 > icat:
                            if cat == cat2:
                                if not cat in excluded:
                                    if not (int(str(catenv[icat].split(":")[1])) == int(
                                            str(catenv[icat2].split(":")[1]))):
                                        hassameCN = False

                                    if not (str(catenv[icat]) == str(catenv[icat2])):
                                        hassameenv = False
                                    if not (connection_corners[icat] == connection_corners[icat2] and connection_edges[
                                        icat] ==
                                            connection_edges[icat2] and connection_faces[icat] == connection_faces[
                                                icat2] and str(catenv[icat]) == str(catenv[icat2])):
                                        # print('here falsch')
                                        hassamechemenvandsameconnectionnumber = False
                                        break

                                else:
                                    print("Test")

        CNsokay = True
        if maxCN != None:
            for icat, cat in enumerate(catid):
                if int(str(catenv[icat].split(":")[1])) > maxCN:
                    CNsokay = False
                    print(catenv[icat])
                    break;
        Outputdict = {}

        # checks if number of cations is equal to the number of unique cations
        notincluded = not (len(uniquecat) != len(catid) and not skip and CNsokay)
        Outputdict['not included in study'] = notincluded
        if not skip:
            Outputdict['hassameCN'] = hassameCN
            Outputdict['hassameenv'] = hassameenv
            Outputdict['hassameenvadsameconnectionnumber'] = hassamechemenvandsameconnectionnumber

        return Outputdict

        # with open('exceptionsCN.json', 'w') as file2:
        #     json.dump(exceptionsCN, file2)
        #
        # with open('exceptionsenv.json', 'w') as file2:
        #     json.dump(exceptionsenv, file2)
        #
        # with open('exceptionsenv_con.json', 'w') as file3:
        #     json.dump(exceptionsenv_con, file3)
        #
        # # TODO maybe good to save names of structures in an extra file the most important informtaion on atoms, coordination environment and number of connections.
        # # also print the information?
        #
        # print('Pauling rule?')
        # print(hassamechemenvandsameconnectionnumber)

        pass


# classes for plotting
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib as mpl

mpl.rcParams["savefig.directory"] = os.chdir(os.getcwd())
mpl.rcParams["savefig.format"] = 'pdf'


class PlotterPSE:
    """

    """

    # was soll das ding koennen
    # plotte z.B. die Ergebnisse in Form des PSE?
    # Mache diesen Plot fuer verschiedene Valenzen bei gleicher Atomsorte und sortiere nach Ordnungszahl

    def __init__(self, cationlist, valuestoplot):
        # valuestoplot is so far a dict with lists for each cation
        # this list contains a positive and negative value
        self.cationlist = cationlist
        self.valuestoplot = valuestoplot

    def _get_fulfilled(self, okay, notokay):
        okayperc = np.float(okay) / np.float(okay + notokay)
        return okayperc

    def set_colormap(self, lowerlimit, upperlimit):
        self.fig.colorbar(self.img1, ax=self.ax).set_clim(lowerlimit, upperlimit)
        # self.ax.clim(lowerlimit, upperlimit)

    def get_plot(self, xlim=[1, 18], ylim=[1, 9]):
        import matplotlib
        # TODO: cleanup
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        matplotlib.rcParams['font.size'] = 3
        mpl.rcParams['figure.dpi'] = 400
        matplotlib.rc("savefig", dpi=400)
        # font = {'family': 'normal',
        #        'size': 10}

        # matplotlib.rc('font', **font)

        PSE = np.full((int(ylim[1] - ylim[0] + 2),
                       int(xlim[1] - xlim[0] + 2)), np.nan)
        # for cat in self.cationlist:
        #     if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
        #         PSE[Element(cat).row, Element(cat).group] = self._get_fulfilled(self.valuestoplot[cat][0],
        #                                                                         self.valuestoplot[cat][1])
        for cat in self.cationlist:
            if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                PSE[Element(cat).row, Element(cat).group] = self._get_fulfilled(self.valuestoplot[cat][0],
                                                                                self.valuestoplot[cat][1])
            else:
                print(cat)
                print(self.valuestoplot[cat][0])
                print(self.valuestoplot[cat][1])
        # print(self.cationlist)
        # neuer plotversuch: jetzt mit farbschema
        fig, ax = plt.subplots()
        # plt.clim(0,1)

        mycm = plt.cm.get_cmap('cool', 24)
        # mycm.set_clim(vmin=0.0, vmax=1.0)
        img1 = ax.imshow(PSE, cmap=mycm)

        # cbar = mpl.colorbar.ColorbarBase(ax, cmap=mycm)
        # cbar.set_clim(0, 1.0)

        # print(mycm.get_clim())

        for cat in self.cationlist:
            if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                ax.text(Element(cat).group - 0.3, Element(cat).row + 0.15, cat)

        ax.set_ylim(np.float(ylim[1]) - 0.5, np.float(ylim[0] - 0.5))
        ax.set_xlim(np.float(xlim[0] - 0.5), np.float(xlim[1] - 0.5))
        ax.set_xticks(range(int(xlim[0]), int(xlim[1]), 1))
        ax.set_yticks(range(int(ylim[0]), int(ylim[1]), 1))
        current_cmap = mpl.cm.get_cmap()
        current_cmap.set_bad(color='white')
        # current_cmap.set_clim(vmin=0.0, vmax=1.0)
        # fig.colorbar(img1, ax=ax).set_clim(0.0, 1.0)
        self.fig = fig
        self.img1 = img1
        self.ax = ax
        # ax.clim(0, 1)
        return plt
        # plt.show()

    def show(self, xlim=[1, 18], ylim=[1, 9]):
        self.get_plot(xlim=xlim, ylim=ylim).show()


class PlotterPSE_Mendeleev(PlotterPSE):
    """

    """
    # was soll das ding koennen
    # plotte z.B. die Ergebnisse in Form des PSE?
    # Mache diesen Plot fuer verschiedene Valenzen bei gleicher Atomsorte und sortiere nach Ordnungszahl

    mendeleev_dict = {
        "H": {"group": 1, "row": 1, "subrow": 1},
        "Li": {"group": 1, "row": 2, "subrow": 2},
        "Na": {"group": 1, "row": 3, "subrow": 3},
        "K": {"group": 1, "row": 4, "subrow": 4},
        "Cu": {"group": 1, "row": 4, "subrow": 5},
        "Rb": {"group": 1, "row": 5, "subrow": 6},
        "Ag": {"group": 1, "row": 5, "subrow": 7},
        "Cs": {"group": 1, "row": 6, "subrow": 8},
        "Au": {"group": 1, "row": 6, "subrow": 9},
        "Fr": {"group": 1, "row": 7, "subrow": 10},
        "Rg": {"group": 1, "row": 7, "subrow": 11},

        "Be": {"group": 2, "row": 2, "subrow": 2},
        "Mg": {"group": 2, "row": 3, "subrow": 3},
        "Ca": {"group": 2, "row": 4, "subrow": 4},
        "Zn": {"group": 2, "row": 4, "subrow": 5},
        "Sr": {"group": 2, "row": 5, "subrow": 6},
        "Cd": {"group": 2, "row": 5, "subrow": 7},
        "Ba": {"group": 2, "row": 6, "subrow": 8},
        "Hg": {"group": 2, "row": 6, "subrow": 9},
        "Ra": {"group": 2, "row": 7, "subrow": 10},
        "Cn": {"group": 2, "row": 7, "subrow": 11},

        "B": {"group": 3, "row": 2, "subrow": 2},
        "Al": {"group": 3, "row": 3, "subrow": 3},
        "Sc": {"group": 3, "row": 4, "subrow": 4},
        "Ga": {"group": 3, "row": 4, "subrow": 5},
        "Y": {"group": 3, "row": 5, "subrow": 6},
        "In": {"group": 3, "row": 5, "subrow": 7},
        "La": {"group": 3, "row": 6, "subrow": 8},
        "Tl": {"group": 3, "row": 6, "subrow": 9},
        "Ac": {"group": 3, "row": 7, "subrow": 10},
        "Nh": {"group": 3, "row": 7, "subrow": 11},

        "C": {"group": 4, "row": 2, "subrow": 2},
        "Si": {"group": 4, "row": 3, "subrow": 3},
        "Ti": {"group": 4, "row": 4, "subrow": 4},
        "Ge": {"group": 4, "row": 4, "subrow": 5},
        "Zr": {"group": 4, "row": 5, "subrow": 6},
        "Sn": {"group": 4, "row": 5, "subrow": 7},
        "Hf": {"group": 4, "row": 6, "subrow": 8},
        "Pb": {"group": 4, "row": 6, "subrow": 9},
        "Rf": {"group": 4, "row": 7, "subrow": 10},
        "Fl": {"group": 4, "row": 7, "subrow": 11},

        "N": {"group": 5, "row": 2, "subrow": 2},
        "P": {"group": 5, "row": 3, "subrow": 3},
        "V": {"group": 5, "row": 4, "subrow": 4},
        "As": {"group": 5, "row": 4, "subrow": 5},
        "Nb": {"group": 5, "row": 5, "subrow": 6},
        "Sb": {"group": 5, "row": 5, "subrow": 7},
        "Ta": {"group": 5, "row": 6, "subrow": 8},
        "Bi": {"group": 5, "row": 6, "subrow": 9},
        "Db": {"group": 5, "row": 7, "subrow": 10},
        "Mc": {"group": 5, "row": 7, "subrow": 11},

        "O": {"group": 6, "row": 2, "subrow": 2},
        "S": {"group": 6, "row": 3, "subrow": 3},
        "Cr": {"group": 6, "row": 4, "subrow": 4},
        "Se": {"group": 6, "row": 4, "subrow": 5},
        "Mo": {"group": 6, "row": 5, "subrow": 6},
        "Te": {"group": 6, "row": 5, "subrow": 7},
        "W": {"group": 6, "row": 6, "subrow": 8},
        "Po": {"group": 6, "row": 6, "subrow": 9},
        "Sg": {"group": 6, "row": 7, "subrow": 10},
        "Lv": {"group": 6, "row": 7, "subrow": 11},

        "F": {"group": 7, "row": 2, "subrow": 2},
        "Cl": {"group": 7, "row": 3, "subrow": 3},
        "Mn": {"group": 7, "row": 4, "subrow": 4},
        "Br": {"group": 7, "row": 4, "subrow": 5},
        "Tc": {"group": 7, "row": 5, "subrow": 6},
        "I": {"group": 7, "row": 5, "subrow": 7},
        "Re": {"group": 7, "row": 6, "subrow": 8},
        "At": {"group": 7, "row": 6, "subrow": 9},
        "Bh": {"group": 7, "row": 7, "subrow": 10},
        "Ts": {"group": 7, "row": 7, "subrow": 11},

        "He": {"group": 8, "row": 1, "subrow": 1},
        "Ne": {"group": 8, "row": 2, "subrow": 2},
        "Ar": {"group": 8, "row": 3, "subrow": 3},
        "Fe": {"group": 8, "row": 4, "subrow": 4},
        "Kr": {"group": 8, "row": 4, "subrow": 5},
        "Ru": {"group": 8, "row": 5, "subrow": 6},
        "Xe": {"group": 8, "row": 5, "subrow": 7},
        "Os": {"group": 8, "row": 6, "subrow": 8},
        "Rn": {"group": 8, "row": 6, "subrow": 9},
        "Hs": {"group": 8, "row": 7, "subrow": 10},
        "Og": {"group": 8, "row": 7, "subrow": 11},

        "Co": {"group": 9, "row": 4, "subrow": 4},
        "Rh": {"group": 9, "row": 5, "subrow": 6},
        "Ir": {"group": 9, "row": 6, "subrow": 8},
        "Mt": {"group": 9, "row": 7, "subrow": 10},

        "Ni": {"group": 10, "row": 4, "subrow": 4},
        "Pd": {"group": 10, "row": 5, "subrow": 6},
        "Pt": {"group": 10, "row": 6, "subrow": 8},
        "Ds": {"group": 10, "row": 7, "subrow": 10},

        "La": {"group": 1, "row": 8, "subrow": 11},
        "Ac": {"group": 1, "row": 9, "subrow": 12},

        "Ce": {"group": 2, "row": 8, "subrow": 11},
        "Th": {"group": 2, "row": 9, "subrow": 12},

        "Pr": {"group": 3, "row": 8, "subrow": 11},
        "Pa": {"group": 3, "row": 9, "subrow": 12},

        "Nd": {"group": 4, "row": 8, "subrow": 11},
        "U": {"group": 4, "row": 9, "subrow": 12},

        "Pm": {"group": 5, "row": 8, "subrow": 11},
        "Np": {"group": 5, "row": 9, "subrow": 12},

        "Sm": {"group": 6, "row": 8, "subrow": 11},
        "Pu": {"group": 6, "row": 9, "subrow": 12},

        "Eu": {"group": 7, "row": 8, "subrow": 11},
        "Am": {"group": 7, "row": 9, "subrow": 12},

        "Gd": {"group": 8, "row": 8, "subrow": 11},
        "Cm": {"group": 8, "row": 9, "subrow": 12},

        "Tb": {"group": 9, "row": 8, "subrow": 11},
        "Bk": {"group": 9, "row": 9, "subrow": 12},

        "Dy": {"group": 10, "row": 8, "subrow": 11},
        "Cf": {"group": 10, "row": 9, "subrow": 12},

        "Ho": {"group": 11, "row": 8, "subrow": 11},
        "Es": {"group": 11, "row": 9, "subrow": 12},

        "Er": {"group": 12, "row": 8, "subrow": 11},
        "Fm": {"group": 12, "row": 9, "subrow": 12},

        "Tm": {"group": 13, "row": 8, "subrow": 11},
        "Md": {"group": 13, "row": 9, "subrow": 12},

        "Yb": {"group": 14, "row": 8, "subrow": 11},
        "No": {"group": 14, "row": 9, "subrow": 12},

        "Lu": {"group": 15, "row": 8, "subrow": 11},
        "Lr": {"group": 15, "row": 9, "subrow": 12},

    }

    def __init__(self, cationlist, valuestoplot):
        # valuestoplot is so far a dict with lists for each cation
        # this list contains a positive and negative value
        self.cationlist = cationlist
        self.valuestoplot = valuestoplot

    def get_mendeleev_group(self, elementstring):
        return self.mendeleev_dict[elementstring]["group"]

    def get_mendeleev_subrow(self, elementstring):
        return self.mendeleev_dict[elementstring]["subrow"]

    def get_mendeleev_mainrow(self, elementstring):
        return self.mendeleev_dict[elementstring]["row"]

    def get_plot(self, xlim=[1, 15], ylim=[1, 12]):
        import matplotlib
        # TODO: cleanup
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        matplotlib.rcParams['font.size'] = 3
        mpl.rcParams['figure.dpi'] = 400
        matplotlib.rc("savefig", dpi=400)
        # font = {'family': 'normal',
        #        'size': 10}

        # matplotlib.rc('font', **font)

        PSE = np.full((int(ylim[1] - ylim[0] + 2),
                       int(xlim[1] - xlim[0] + 2)), np.nan)
        # for cat in self.cationlist:
        #     if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
        #         PSE[Element(cat).row, Element(cat).group] = self._get_fulfilled(self.valuestoplot[cat][0],
        #                                                                         self.valuestoplot[cat][1])
        for cat in self.cationlist:
            if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                print(cat)
                PSE[self.get_mendeleev_subrow(str(cat)), self.get_mendeleev_group(str(cat))] = self._get_fulfilled(
                    self.valuestoplot[cat][0],
                    self.valuestoplot[cat][1])
            else:
                print(cat)
                print(self.valuestoplot[cat][0])
                print(self.valuestoplot[cat][1])
        # print(self.cationlist)
        # neuer plotversuch: jetzt mit farbschema
        fig, ax = plt.subplots()
        # plt.clim(0,1)

        mycm = plt.cm.get_cmap('cool', 24)
        # mycm.set_clim(vmin=0.0, vmax=1.0)
        img1 = ax.imshow(PSE, cmap=mycm)

        # cbar = mpl.colorbar.ColorbarBase(ax, cmap=mycm)
        # cbar.set_clim(0, 1.0)

        # print(mycm.get_clim())

        for cat in self.cationlist:
            if not (self.valuestoplot[cat][0] == 0 and self.valuestoplot[cat][1] == 0):
                ax.text(self.get_mendeleev_group(str(cat)) - 0.3, self.get_mendeleev_subrow(str(cat)) + 0.15, cat)

        ax.set_ylim(np.float(ylim[1]) - 0.5, np.float(ylim[0] - 0.5))
        ax.set_xlim(np.float(xlim[0] - 0.5), np.float(xlim[1] - 0.5))
        ax.set_xticks(range(int(xlim[0]), int(xlim[1]), 1))
        ax.set_yticks(range(int(ylim[0]), int(ylim[1]), 1))
        current_cmap = mpl.cm.get_cmap()
        current_cmap.set_bad(color='white')
        # current_cmap.set_clim(vmin=0.0, vmax=1.0)
        # fig.colorbar(img1, ax=ax).set_clim(0.0, 1.0)
        self.fig = fig
        self.img1 = img1
        self.ax = ax
        # ax.clim(0, 1)
        return plt
        # plt.show()

    def show(self, xlim=[1, 18], ylim=[1, 9]):
        self.get_plot(xlim=xlim, ylim=ylim).show()


class PlotterZdependent:
    def __init__(self, cationlist, cat_valence_list, sort_by_Z=True, remove_cat_only_one_val=True, lowerlimit=0):
        # sort cat_valence_list first by Z

        if sort_by_Z:
            sortedcationlist = []
            for el in Element:
                if str(el) in cationlist:
                    sortedcationlist.append(str(el))
            cationlist = sortedcationlist
            pass
            # sortedcationlist=[]
            # for Z in range(1,103):
            #     #if str(Element(Z)) in cationlist:
            #     #    sortedcationlist.append(str(Element(Z)))

        self.catswithval = []
        # self.SIZE=len(cat_valence_list)
        MAXVALENCE = 0
        for icat, cat in enumerate(cationlist):
            localmax = 0
            for key in cat_valence_list[cat].keys():
                if localmax < int(key):
                    localmax = int(key)
            if localmax > MAXVALENCE:
                MAXVALENCE = localmax

        self.MAXVALENCE = MAXVALENCE

        for icat, cat in enumerate(cationlist):
            hasmorethanonevalence = False
            countvalences = 0
            for val in range(0, MAXVALENCE + 1):
                if val in cat_valence_list[cat]:
                    if cat_valence_list[cat][val] != [0, 0]:
                        countvalences = countvalences + 1

            if (not (countvalences == 1 or countvalences == 0)) and remove_cat_only_one_val:
                self.catswithval.append(cat)
            elif (countvalences != 0) and (not remove_cat_only_one_val):
                self.catswithval.append(cat)
        self.valenzabh = np.full((MAXVALENCE + 1, len(self.catswithval)), np.nan)
        for val in range(0, MAXVALENCE + 1):
            for icat, cat in enumerate(self.catswithval):
                if val in cat_valence_list[cat]:
                    if cat_valence_list[cat][val] != [0, 0]:
                        if (cat_valence_list[cat][val][0] + cat_valence_list[cat][val][1]) > lowerlimit:
                            self.valenzabh[val, icat] = self._get_fulfilled(cat_valence_list[cat][val][0],
                                                                            cat_valence_list[cat][val][1])
        self.SIZE = len(self.catswithval)

    def _get_fulfilled(self, okay, notokay):
        okayperc = okay / np.float(okay + notokay)
        return okayperc

    def get_plot(self):
        import matplotlib
        matplotlib.rcParams['pdf.fonttype'] = 42
        matplotlib.rcParams['ps.fonttype'] = 42
        matplotlib.rcParams['font.size'] = 2

        fig, ax = plt.subplots()
        img1 = ax.imshow(self.valenzabh, origin="lower", cmap=plt.cm.get_cmap('cool', 24))
        ax.set_xticks(range(0, self.SIZE + 1))
        ax.set_xlim(-0.5, self.SIZE - 0.5)
        ax.set_xticklabels(self.catswithval)
        ax.set_ylim(0.5, self.MAXVALENCE + 0.5)
        ax.set_yticks(range(1, self.MAXVALENCE + 1))
        ax.set_ylabel('Valence')
        current_cmap = mpl.cm.get_cmap()
        current_cmap.set_bad(color='white')
        fig.colorbar(img1, ax=ax)
        return plt

    def show(self):
        self.get_plot()
        plt.show()


def get_compounds_lists():
    ce_fraction = 0.95
    csm = 2
    anion = 'O'
    charge = -2
    e_above_hull = 0.025
    filename = 'ce_fraction_' + \
               str(ce_fraction) + 'plus_csm_' + str(csm) + 'plus_eabovehull' + str(e_above_hull) + '.json'

    # to get it from the database:
    Materials = is_clear_materials()
    # list_compound=Materials.get_isclearmaterials_fromdatabase(ce_fraction=ce_fraction, csm=csm,anion=anion,charge=charge,e_above_hull=e_above_hull)
    # Materials.save_infile(ce_fraction=ce_fraction,csm=csm,anion=anion,charge=charge,filename=filename,e_above_hull=e_above_hull)

    list_compound = Materials.load_is_clear_materials_fromfile(filename=filename, ce_fraction=ce_fraction, csm=csm,
                                                               anion=anion, charge=charge, e_above_hull=e_above_hull)

    print(list_compound)

    list_compound2 = get_materials_standard()

    list_compound.extend(list_compound2)

    print(len(set(list_compound)))


def main():
    #    pass
    # secondrule_analyze_more_detailed()

    # fourthrule_plotresults_CN_valence_sep()
    # thirdrule_structureanalysis_plot()
    # get_compounds_lists()

    # thirdrule_plotresults_CN_valence_sep()

    # extend this ratio stuff
    # combine with and without shannon radii

    firstrule(ratiosaccordingtoPauling=True, onlylowerlimit=True, e_above_hull=0.025)
# firstrule_withShannon(ratiosaccordingtoPauling=False,onlylowerlimit=False)
# plots have been tested

# TODO: make prediction with mean Shannon radii?

# TODO: Einheitliche Zahl Polyeder fuer erste Regel verwendet?
# so far there is not lower limit for the number of cations considered
# Pruefen
# firstrulecheckup(ratiosaccordingtoPauling=True)
# secondrule()
# TODO: Diese hier genau testen und Plots anschauen
# Auch Tortenplots machen
# thirdrule()
# thirdrule_analysestruct()
# TODO: Skala von 0.0 bis 1.0 ausdehnen
#
# diesen teil in andere third rule einbeziehen
# thirdrule_plotresults_each()
## thirdrule_plotresults_worked()
# fourthrule()
# fourthrule_focusonpolyhedra()
# fourthrule_plotresults_worked2()

# fifthrule()

# fifthruleaddition()

# main()

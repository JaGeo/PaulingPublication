import json
import numpy as np
from pymatgen.core import PeriodicSite
import os


class RuleCannotBeAnalyzedError(Exception):
    def __init__(self, value='The Rule cannot be analyzed'):
        self.value = value

    def __str__(self):
        return repr(self.value)


def is_an_oxide_and_no_env_for_O(lse):
    for isite, site in enumerate(lse.structure):
        if lse.valences[isite] < 0 and site.species_string != 'O':
            raise ValueError("This is not an oxide. The assessment will be stopped.")
    for isite, site_envs in enumerate(lse.coordination_environments):
        if lse.structure[isite].species_string=='O':
            if site_envs!=None:
                raise ValueError("Site_envs of anions have been computed. The code has to stop. Use only_cations in compute_structure_environments")
    return True


class Pauling1:
    def __init__(self, lse, filenameradii="univalent_cat_radii.json", onlylowerlimit=False):
        """
        class to test Pauling's first rule
        :param lse: LightStructureEnvironment
        :param filenameradii: name of the file containing radius ratios
        :param onlylowerlimit: If False, ratio windows are considered
        """
        is_an_oxide_and_no_env_for_O(lse)

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
            raise RuleCannotBeAnalyzedError("The first rule cannot be evaluated.")

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
        is_an_oxide_and_no_env_for_O(lse)

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

    def newsetup(self, lse, filename=None, save_to_file=True, foldername='ThirdRuleAnalysisConnections', distance=8.0):
        """
            :param lse: LightStructureEnvironment
            :param save_to_file: Boolean
            :param filename: beginning of the file name, without ".json"
            :param foldername: name of the folder
            :param distance: float giving the distances of cations that is considered
        """
        is_an_oxide_and_no_env_for_O(lse)
        super().__init__(DISTANCE=distance)
        if filename is None:
            save_to_file == False
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
                with open(os.path.join(foldername, self.mat + '.json'), 'w') as file:
                    json.dump(self.PolyhedronDict, file)

    def fromfile(self, filename, foldername='ThirdRuleAnalysisConnections'):
        #TODO: make the use of filename consistent!
        #make sure you can use it like this
        # write it differently
        self.mat = filename
        with open(foldername + '/' + mat + '.json') as file:
            self.PolyhedronDict = json.load(file)

    def _test_fulfillment(self, maxCN=None):
        """
        internal class to test fulfillment of rule
        :param maxCN:
        :return: returns connections as a dict
        """
        inputdict = self.PolyhedronDict
        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        not_connected = 0
        corner = 0
        edge = 0
        face = 0
        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if maxCN is None:
                if herepolyhedra[numberpolyhedra] == 0:
                    not_connected += 1
                if herepolyhedra[numberpolyhedra] == 1:
                    corner += 1
                if herepolyhedra[numberpolyhedra] == 2:
                    edge += 1
                if herepolyhedra[numberpolyhedra] >= 3:
                    face += 1
            else:
                if (iinfo['CN'][0] <= maxCN) and (iinfo['CN'][1] <= maxCN):
                    if herepolyhedra[numberpolyhedra] == 0:
                        not_connected += 1
                    if herepolyhedra[numberpolyhedra] == 1:
                        corner += 1
                    if herepolyhedra[numberpolyhedra] == 2:
                        edge += 1
                    if herepolyhedra[numberpolyhedra] >= 3:
                        face += 1

            numberpolyhedra += 1

        return {"no": not_connected, "corner": corner, "edge": edge, "face": face}

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
                # inspired by Lattice class in pymatgen.core.lattice
                recp_len = np.array(struct.lattice.reciprocal_lattice.abc) / (2 * np.pi)
                nmax = float(DISTANCE) * recp_len + 0.01

                # +1 due to the fact that also neighbors of an atom on (1,1,1) have to be considered
                ceila = int(np.ceil(nmax[0]) + 1)
                ceilb = int(np.ceil(nmax[1]) + 1)
                ceilc = int(np.ceil(nmax[2]) + 1)

                for i in range(0, ceila):
                    for j in range(0, ceilb):
                        for k in range(0, ceilc):
                            addcoord = [float(i), float(j), float(k)]
                            newcoords = coords + addcoord
                            allsites.append(PeriodicSite(
                                atoms_n_occu, newcoords, lattice, to_unit_cell=False))

        return allsites

    def get_connections(self, maximumCN=None):
        """
        Gives connections
        :return: Outputdict with number of connections
        """
        OutputDict = self._test_fulfillment(maxCN=maximumCN)
        return OutputDict


class Pauling3(Pauling3and4):
    def is_fulfilled(self, maximumCN=None):
        """
        tells you if third rule is fulfilled (i.e. no connections via faces!)
        :param maximumCN: gives the maximal CN of cations that are considered
        :return:
        """
        connections = self._test_fulfillment(maxCN=maximumCN)

        if connections["face"] == 0:
            return True
        else:
            return False

    def get_details(self, maximumCN=None):
        """
        returns a dict that looks the following {"Element string": "Valence as an int": {"no": number1, "corner": number2, "edge": number3, "face" number4} indicating the connections.}
        :param maximumCN: maximum CN that is considered
        :return: Outputdict with information on the types of connetions of pairs of polyhedra and with information on the valences, species that are connected via corners, edges, and faces or not connected

        """
        Outputdict = self._postevaluation3rdrule_evaluate_elementdependency(maxCN=maximumCN)
        seconddict = self._test_fulfillment(maxCN=maximumCN)
        seconddict["species"] = Outputdict

        return seconddict

    def _postevaluation3rdrule_evaluate_elementdependency(self, maxCN=None):
        """
        evaluates element and valence dependency of connected pairs of polyhedra
        :param maxCN:
        :return: a Outputdict that will be evaluated by other methods
        """
        inputdict = self.PolyhedronDict
        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        # Dict that includes cations sharing no connections, corners, edges, and faces
        # This can be used to evaluate the results per polyhedron pair
        Dict_ElementDependency = {}
        numberpolyhedra = 0
        for iinfo in additionalinfo:
            if maxCN is None:
                if iinfo['cations'][0] not in Dict_ElementDependency:
                    Dict_ElementDependency[iinfo['cations'][0]] = {}
                if iinfo['cations'][1] not in Dict_ElementDependency:
                    Dict_ElementDependency[iinfo['cations'][1]] = {}

                if iinfo['valences'][0] not in Dict_ElementDependency[iinfo['cations'][0]]:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]] = {'no': 0, 'corner': 0,
                                                                                         'edge': 0, 'face': 0}

                if iinfo['valences'][1] not in Dict_ElementDependency[iinfo['cations'][1]]:
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]] = {'no': 0, 'corner': 0,
                                                                                         'edge': 0, 'face': 0}

                if herepolyhedra[numberpolyhedra] == 0:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['no'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['no'] += 1

                elif herepolyhedra[numberpolyhedra] == 1:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['corner'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['corner'] += 1


                elif herepolyhedra[numberpolyhedra] == 2:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['edge'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['edge'] += 1


                elif herepolyhedra[numberpolyhedra] > 2:
                    Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['face'] += 1
                    Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['face'] += 1

            else:
                if iinfo["CN"][0] <= maxCN and iinfo["CN"][1] <= maxCN:
                    if iinfo['cations'][0] not in Dict_ElementDependency:
                        Dict_ElementDependency[iinfo['cations'][0]] = {}
                    if iinfo['cations'][1] not in Dict_ElementDependency:
                        Dict_ElementDependency[iinfo['cations'][1]] = {}

                    if iinfo['valences'][0] not in Dict_ElementDependency[iinfo['cations'][0]]:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]] = {'no': 0, 'corner': 0,
                                                                                             'edge': 0, 'face': 0}

                    if iinfo['valences'][1] not in Dict_ElementDependency[iinfo['cations'][1]]:
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]] = {'no': 0, 'corner': 0,
                                                                                             'edge': 0, 'face': 0}

                    if herepolyhedra[numberpolyhedra] == 0:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['no'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['no'] += 1

                    elif herepolyhedra[numberpolyhedra] == 1:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['corner'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['corner'] += 1


                    elif herepolyhedra[numberpolyhedra] == 2:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['edge'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['edge'] += 1


                    elif herepolyhedra[numberpolyhedra] > 2:
                        Dict_ElementDependency[iinfo['cations'][0]][iinfo['valences'][0]]['face'] += 1
                        Dict_ElementDependency[iinfo['cations'][1]][iinfo['valences'][1]]['face'] += 1

            numberpolyhedra = numberpolyhedra + 1

        return Dict_ElementDependency


class Pauling4(Pauling3and4):
    def is_fulfilled(self):
        """
        tells you if polyhedra of cations with highest valence and smallest CN don't show any connections within the structure
        structure has to have cations with different valences and coordination numbers
        :return: Boolean
        :raises: RuleCannotBeAnalyzedError if structure cannot will be accessed for this test (no cation that differ in valence and CN)
        """
        if not self._is_candidate_4thrule():
            raise RuleCannotBeAnalyzedError("The fourth rule cannot be evaluated")

        maxval = max(self.PolyhedronDict['cationvalences'])
        minCN = min(self.PolyhedronDict['CNlist'])
        outputdict = self._postevaluation4thruleperpolyhedron_only_withoutproduct(minCN, minCN, maxval, maxval)
        if outputdict['corner'] != 0 or outputdict['edge'] != 0 or outputdict['face'] != 0:
            return False
        else:
            return True

    def get_details(self):
        """
        gives you number of connections as a function of the coordination numbers and valences of the polyhedra pairs
        :return: Dict of the following form: {"val1:valence1": {"val2:valence2": {"CN1:CN1": {"CN2:CN2": {"no": number1, "corner": number2, "edge": number3, "face" number4}}}}} indicating the connections depending on valences and CN
        """
        if not self._is_candidate_4thrule():
            raise RuleCannotBeAnalyzedError()
        return self._postevaluation4thrule()

    def _is_candidate_4thrule(self):
        inputdict = self.PolyhedronDict
        samevalences = inputdict['samevalences']
        sameCN = inputdict['sameCN']
        if (not samevalences) or (not sameCN):
            return True
        else:
            return False

    # Add additional evaluation of 4th rule that does not depend on valence or CN product but on 2 valences, 2 CN
    def _postevaluation4thrule(self):
        """
        :return: more information connected pairs of polyhedra
        """
        inputdict = self.PolyhedronDict

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']

        Outputdict = {}
        Elementwise = {}
        numberpolyhedra = 0
        for iinfo in additionalinfo:

            if not iinfo['cations'][0] in Elementwise:
                Elementwise[iinfo['cations'][0]] = {}
            if not iinfo['cations'][1] in Elementwise:
                Elementwise[iinfo['cations'][1]] = {}

            # here: consider elements:
            # symmetry of valence and CN have to be considered
            if not ("val1:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][0]]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])] = {}
            if not ("val1:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][0]]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][0]][
                "val1:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][0]][
                "val1:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])] = {}

            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])] = {}

            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}

            if herepolyhedra[numberpolyhedra] == 0:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
            elif herepolyhedra[numberpolyhedra] == 1:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1

            elif herepolyhedra[numberpolyhedra] == 2:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1

            elif herepolyhedra[numberpolyhedra] > 2:
                Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (info['valences'][0] == info['valences'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (info['CN'][0] == info['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1
                if not (info['valences'][0] == info['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][0]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1

            # here: consider elements:
            # symmetry of valence and CN have to be considered
            if not ("val1:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][1]]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])] = {}
            if not ("val1:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][1]]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][1])) in Elementwise[iinfo['cations'][1]][
                "val1:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][0])) in Elementwise[iinfo['cations'][1]][
                "val1:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])] = {}

            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])] = {}

            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                       "val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                       "val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                    "val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}

            if herepolyhedra[numberpolyhedra] == 0:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
            elif herepolyhedra[numberpolyhedra] == 1:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1

            elif herepolyhedra[numberpolyhedra] == 2:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1

            elif herepolyhedra[numberpolyhedra] > 2:
                Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                    "val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (info['valences'][0] == info['valences'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (info['CN'][0] == info['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][0])][
                        "val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1
                if not (info['valences'][0] == info['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Elementwise[iinfo['cations'][1]]["val1:" + str(iinfo['valences'][1])][
                        "val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1

            # symmetry of valence and CN have to be considered
            if not ("val1:" + str(iinfo['valences'][0])) in Outputdict:
                Outputdict["val1:" + str(iinfo['valences'][0])] = {}
            if not ("val1:" + str(iinfo['valences'][1])) in Outputdict:
                Outputdict["val1:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][1])) in Outputdict["val1:" + str(iinfo['valences'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])] = {}
            if not ("val2:" + str(iinfo['valences'][0])) in Outputdict["val1:" + str(iinfo['valences'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])] = {}

            if not ("CN1:" + str(iinfo['CN'][0])) in Outputdict["val1:" + str(iinfo['valences'][0])][
                "val2:" + str(iinfo['valences'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][0])) in Outputdict["val1:" + str(iinfo['valences'][1])][
                "val2:" + str(iinfo['valences'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in Outputdict["val1:" + str(iinfo['valences'][0])][
                "val2:" + str(iinfo['valences'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])] = {}
            if not ("CN1:" + str(iinfo['CN'][1])) in Outputdict["val1:" + str(iinfo['valences'][1])][
                "val2:" + str(iinfo['valences'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])] = {}

            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][1])) in \
                   Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][0])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}
            if not ("CN2:" + str(iinfo['CN'][0])) in \
                   Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                       "CN1:" + str(iinfo['CN'][1])]:
                Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                    "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])] = {"no": 0, "corner": 0, "edge": 0,
                                                                                   "face": 0}

            if herepolyhedra[numberpolyhedra] == 0:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["no"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["no"] += 1
            elif herepolyhedra[numberpolyhedra] == 1:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["corner"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["corner"] += 1

            elif herepolyhedra[numberpolyhedra] == 2:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["edge"] += 1
                if not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1
                if not (iinfo['valences'][0] == iinfo['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["edge"] += 1

            elif herepolyhedra[numberpolyhedra] > 2:
                Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                    "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (info['valences'][0] == info['valences'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][0])]["CN2:" + str(iinfo['CN'][1])]["face"] += 1
                if not (info['CN'][0] == info['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][0])]["val2:" + str(iinfo['valences'][1])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1
                if not (info['valences'][0] == info['valences'][1]) or not (iinfo['CN'][0] == iinfo['CN'][1]):
                    Outputdict["val1:" + str(iinfo['valences'][1])]["val2:" + str(iinfo['valences'][0])][
                        "CN1:" + str(iinfo['CN'][1])]["CN2:" + str(iinfo['CN'][0])]["face"] += 1

            numberpolyhedra = numberpolyhedra + 1

        Outputdict['elementwise'] = Elementwise

        return Outputdict

    def _postevaluation4thruleperpolyhedron_only_withoutproduct(self, CN1, CN2, val1, val2):
        """
        That is the one I used for the evaluation
        :param CN1: CN1 of an atom in pair of polyhedra
        :param CN2: CN2 of the other atom in pair of polyhedra
        :param val1: valence 1 of an atom in pair of polyhedra
        :param val2: valence 2 of an atom in pair of polyhedra
        :return: more information connected pairs of polyhedra
        """
        inputdict = self.PolyhedronDict

        additionalinfo = inputdict['Additional']
        herepolyhedra = inputdict['PolyConnect']
        samevalences = inputdict['samevalences']
        sameCN = inputdict['sameCN']

        notconnected = 0
        corner = 0
        edge = 0
        face = 0

        numberpolyhedra = 0
        for iinfo in additionalinfo:
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
        Thirdruledict['no'] = notconnected
        Thirdruledict['corner'] = corner
        Thirdruledict['edge'] = edge
        Thirdruledict['face'] = face

        return Thirdruledict


class Pauling5(PaulingConnection):

    def __init__(self):
        pass

    def newsetup(self, lse, filename=None, save_to_file=True, foldername="FifthsRuleAnalysis", distance=8.0):

        # collects the coordination numbers, coordination environments and the number of connections via corners, edges and faces of each of the polyhedra for each of the cations
        """
                    :param lse: LightStructureEnvironment
                    :param save_to_file: Boolean
                    :param filename: beginning of the file name, without ".json"
                    :param foldername: name of the folder
                    :param distance: float giving the distances of cations that is considered
        """
        is_an_oxide_and_no_env_for_O(lse)
        super().__init__(DISTANCE=distance)

        if filename is None:
            save_to_file == False
        if save_to_file:
            if not os.path.isdir(foldername):
                os.mkdir(foldername)
        self.mat = filename

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

        # unique cations with valences
        uniquecat = []
        for cat in catid:
            if cat not in uniquecat:
                uniquecat.append(cat)

        outputdict = {}
        outputdict['connection_corners'] = connection_corners
        outputdict['connection_edges'] = connection_edges
        outputdict['connection_faces'] = connection_faces
        outputdict['catenv'] = catenv
        outputdict['catid'] = catid
        outputdict['uniquecat'] = uniquecat

        if save_to_file:
            with open(os.path.join(foldername, self.mat + '.json'), 'w') as file:
                json.dump(outputdict, file)

        self.FifthRuleDict = outputdict

    def fromfile(self, mat, foldername):
        # have to test this part of the code as well
        self.mat = mat
        print(mat)
        with open(foldername + '/' + mat + '.json') as file:
            self.FifthRuleDict = json.load(file)

    def is_fulfilled(self, options="CN"):
        """
        tests if all chemically equivalent cations (same element, same valence) have the same CN, or environment, or environment and number of connections
        if there is only one chemically equivalent cation, the method raises a RuleCannotBeAnalyzedError
        if the wrong option is used, a ValueError is raised
        :param options: can be "CN", "env", or "env+nconnections"
        :return: Boolean
        """
        # test _is_candidate_5thrule
        if not self._is_candidate_5thrule():
            raise RuleCannotBeAnalyzedError("5th Rule cannot be evaluated")

        outputdict = self._postevaluation5thrule()
        if options == 'CN':
            return outputdict['hassameCN']
        elif options == 'env':
            return outputdict['hassameenv']
        elif options == 'env+nconnections':
            return outputdict['hassameenvadsameconnectionnumber']
        else:
            raise ValueError("Wrong Option")

    def get_details(self, options='CN'):
        if not self._is_candidate_5thrule():
            raise RuleCannotBeAnalyzedError("5th Rule cannot be evaluated")

        outputdict = self._postevaluation5thrule_elementdependency()
        output = {}
        if options == 'CN':
            for cat in outputdict['exceptionsCN']:
                if not cat[0] in output:
                    output[cat[0]] = {}
                    output[cat[0]]['not_fulfilled'] = 0
                    output[cat[0]]['fulfilled'] = 0
                output[cat[0]]['not_fulfilled'] += 1
            for cat in outputdict['fulfillingCN']:
                if not cat[0] in output:
                    output[cat[0]] = {}
                    output[cat[0]]['not_fulfilled'] = 0
                    output[cat[0]]['fulfilled'] = 0
                output[cat[0]]['fulfilled'] += 1
            return output
        elif options == 'env':
            for cat in outputdict['exceptionsenvs']:
                if not cat[0] in output:
                    output[cat[0]] = {}
                    output[cat[0]]['not_fulfilled'] = 0
                    output[cat[0]]['fulfilled'] = 0
                output[cat[0]]['not_fulfilled'] += 1
            for cat in outputdict['fulfillingenvs']:
                if not cat[0] in output:
                    output[cat[0]]['not_fulfilled'] = 0
                    output[cat[0]]['fulfilled'] = 0
                output[cat[0]]['fulfilled'] += 1
            return output
        elif options == 'env+nconnections':
            for cat in outputdict['exceptionsenvs_connections']:
                if not cat[0] in output:
                    output[cat[0]] = {}
                    output[cat[0]]['not_fulfilled'] = 0
                    output[cat[0]]['fulfilled'] = 0
                output[cat[0]]['not_fulfilled'] += 1
            for cat in outputdict['fulfillingenvs_connections']:
                if not cat[0] in output:
                    output[cat[0]] = {}
                    output[cat[0]]['not_fulfilled'] = 0
                    output[cat[0]]['fulfilled'] = 0
                output[cat[0]]['fulfilled'] += 1
            return output
        else:
            raise ValueError("Wrong option")

    def _is_candidate_5thrule(self):
        """
        tells you if 5th rule can be evaluated
        :return: Boolean
        """
        catid = self.FifthRuleDict['catid']
        uniquecat = self.FifthRuleDict['uniquecat']
        if (len(uniquecat) != len(catid)):
            return True
        else:
            return False

    def _postevaluation5thrule(self):
        """
        evaluates the 5th rule
        :return: outputdict with relevant information
        """
        connection_corners = self.FifthRuleDict['connection_corners']
        connection_edges = self.FifthRuleDict['connection_edges']
        connection_faces = self.FifthRuleDict['connection_faces']
        catenv = self.FifthRuleDict['catenv']
        catid = self.FifthRuleDict['catid']
        uniquecat = self.FifthRuleDict['uniquecat']

        hassamechemenvandsameconnectionnumber = True
        hassameenv = True
        hassameCN = True

        # tells you if all cations that are the same have the same CN, CE, CE with no connections
        # is this the correct test?
        if len(uniquecat) != len(catid):

            for icat, cat in enumerate(catid):
                for icat2, cat2 in enumerate(catid):
                    if icat2 > icat:
                        if cat == cat2:
                            if not (int(str(catenv[icat].split(":")[1])) == int(
                                    str(catenv[icat2].split(":")[1]))):
                                hassameCN = False

                            if not (str(catenv[icat]) == str(catenv[icat2])):
                                hassameenv = False
                            if not (connection_corners[icat] == connection_corners[icat2] and connection_edges[
                                icat] ==
                                    connection_edges[icat2] and connection_faces[icat] == connection_faces[
                                        icat2] and str(catenv[icat]) == str(catenv[icat2])):
                                hassamechemenvandsameconnectionnumber = False
                                break

        Outputdict = {}
        Outputdict['hassameCN'] = hassameCN
        Outputdict['hassameenv'] = hassameenv
        Outputdict['hassameenvadsameconnectionnumber'] = hassamechemenvandsameconnectionnumber

        return Outputdict

    def _postevaluation5thrule_elementdependency(self):
        """
        gives you information on the element dependency of the rule
        :return:
        """

        connection_corners = self.FifthRuleDict['connection_corners']
        connection_edges = self.FifthRuleDict['connection_edges']
        connection_faces = self.FifthRuleDict['connection_faces']
        catenv = self.FifthRuleDict['catenv']
        catid = self.FifthRuleDict['catid']
        uniquecat = self.FifthRuleDict['uniquecat']

        hassamechemenvandsameconnectionnumber = True
        hassameenv = True
        hassameCN = True

        exceptions_CN = []
        exceptions_envs = []
        exceptions_envs_connections = []

        if len(uniquecat) != len(catid):

            for icat, cat in enumerate(catid):
                for icat2, cat2 in enumerate(catid):
                    if icat2 > icat:
                        if cat == cat2:
                            if not (int(str(catenv[icat].split(":")[1])) == int(
                                    str(catenv[icat2].split(":")[1]))):
                                hassameCN = False
                                if cat not in exceptions_CN:
                                    exceptions_CN.append(cat)
                            if not (str(catenv[icat]) == str(catenv[icat2])):
                                hassameenv = False
                                if cat not in exceptions_envs:
                                    exceptions_envs.append(cat)
                            if not (connection_corners[icat] == connection_corners[icat2] and connection_edges[
                                icat] ==
                                    connection_edges[icat2] and connection_faces[icat] == connection_faces[
                                        icat2] and str(catenv[icat]) == str(catenv[icat2])):
                                # print('here falsch')
                                hassamechemenvandsameconnectionnumber = False

                                if cat not in exceptions_envs_connections:
                                    exceptions_envs_connections.append(cat)

        fulfilling_CN = []
        fulfilling_envs = []
        fulfilling_envs_connections = []
        for unicat in uniquecat:
            if unicat not in exceptions_CN:
                fulfilling_CN.append(unicat)
            if unicat not in exceptions_envs:
                fulfilling_envs.append(unicat)
            if unicat not in exceptions_envs_connections:
                fulfilling_envs_connections.append(unicat)

        Outputdict = {}

        # checks if number of cations is equal to the number of unique cations
        # find a better way for this CNsokay part
        Outputdict['hassameCN'] = hassameCN
        Outputdict['hassameenv'] = hassameenv
        Outputdict['hassameenvadsameconnectionnumber'] = hassamechemenvandsameconnectionnumber
        Outputdict['exceptionsCN'] = exceptions_CN
        Outputdict['exceptionsenvs'] = exceptions_envs
        Outputdict['exceptionsenvs_connections'] = exceptions_envs_connections
        Outputdict['fulfillingCN'] = fulfilling_CN
        Outputdict['fulfillingenvs'] = fulfilling_envs
        Outputdict['fulfillingenvs_connections'] = fulfilling_envs_connections

        return Outputdict

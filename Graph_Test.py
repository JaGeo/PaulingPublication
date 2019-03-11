from PaulingRules import PaulingConnection, is_an_oxide_and_no_env_for_O

class Graph_Connection(PaulingConnection):

    def __init__(self):
        pass

    def newsetup(self, lse, distance=8.0):

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


        struct = lse.structure
        sites = struct.sites
        catid = []
        catenv = []
        connection_corners = []
        connection_edges = []
        connection_faces = []
        Result_Dict={}
        for isite, site in enumerate(sites):
            if self._is_cationic_site(isite, valences=lse.valences):
                Result_Dict[isite] = {'corner': [], 'edge': [], 'face': []}
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
                        Result_Dict[isite]['corner'].append(indexneigh)
                    if numberconnect == 2:
                        edge = edge + 1
                        Result_Dict[isite]['edge'].append(indexneigh)
                    if numberconnect >= 3:
                        face = face + 1
                        Result_Dict[isite]['face'].append(indexneigh)

                connection_corners.append(corner)
                connection_edges.append(edge)
                connection_faces.append(face)
        print(Result_Dict)

        #test which isites have same neighbors
        for isite,connections in Result_Dict.items():
            pass
            #get species for each site and enviroment, number of connections to which kind of elements and sort by this

        #then, check second sphere and so on
        #go to second shell



        #at the end: check if you only find the kind of connections described in the graph for each site

This data was generated with pymatgen version v2019.12.3.
valences = #valences from Bondvalence analyser or the ICSD

#rest of the code!
lgf = LocalGeometryFinder()
lgf.setup_structure(structure=struct)
se = lgf.compute_structure_environments(only_cations=True, valences=valences)
strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()
lse = LightStructureEnvironments.from_structure_environments(strategy=strategy, structure_environments=se)

The structures can be identified with the Materials Project ID. The data is licensed under a CC BY 4.0 license.

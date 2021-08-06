import asterid as ad
def asterid_dm_to_dendropy_dm(D, ts):
    pdm = dendropy.PhylogeneticDistanceMatrix()
    pdm.taxon_namespace = dendropy.TaxonNamespace()
    pdm._mapped_taxa = set()

    for i in range(len(ts)):
        for j in enumerate(ts):
            si = ts[i]
            sj = ts[j]
            dij = D[i, j]

            xi = pdm.taxon_namespace.get_taxon(si)
            if not xi:
                xi = dendropy.Taxon(si)
                pdm.taxon_namespace.add_taxon(xi)
                pdm._mapped_taxa.add(xi)
                pdm._taxon_phylogenetic_distances[xi] = {}

            xj = pdm.taxon_namespace.get_taxon(sj)
            if not xj:
                xj = dendropy.Taxon(sj)
                pdm.taxon_namespace.add_taxon(xj)
                pdm._mapped_taxa.add(xj)
                pdm._taxon_phylogenetic_distances[xj] = {}

            dij = float(dij)
            pdm._taxon_phylogenetic_distances[xi][xj] = dij
    return pdm


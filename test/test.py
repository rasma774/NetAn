import NetAn
import numpy as np

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2022 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'


TEST_GENES = ['IRS1', 'INS', 'INSR']

def test_input():
    nt = NetAn.NetAn(TEST_GENES)
    assert np.all(nt.genes == TEST_GENES)
    assert nt.genes_all is None
    assert nt.genes_exp is None

def test_expansion():
    nt = NetAn.NetAn(TEST_GENES)
    nt.add_close_genes(augment_factor=1)
    assert len(nt.genes) == len(nt.genes_exp)

    nt.add_close_genes(augment_factor=2)
    assert len(nt.genes)*2 == len(nt.genes_exp)

    assert len(nt.genes_exp) == len(np.unique(nt.genes_exp))

def test_enrichments():
    nt = NetAn.NetAn(TEST_GENES)
    nt.add_close_genes(augment_factor=10)

    enr = nt.enrich(gene_set='all', n_clusters=2)
    assert enr['cluster_0'][0].OR.max() > 1

    enr = nt.enrich(gene_set='all', connected_components=True)
    assert enr['component_0'][0].OR.max() > 1

    enr = nt.enrich(gene_set='all', use_cliques=True)
    assert enr[list(enr.keys())[0]][0].OR.max() > 1


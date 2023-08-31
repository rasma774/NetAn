import NetAn
import pandas as pd

MS_genes = pd.read_csv('MS_genes.csv').values.T[0]

nw_annot = NetAn.NetAn(MS_genes, nw_cutoff=.5)

MS_cluster_exp = nw_annot.enrich(gene_set='start_genes',
                                 n_clusters=2,
                                 )

f_exp, ax_exp = NetAn.plot_utils.plot_enrichment_clusters(MS_cluster_exp, sorton='pval')

MS_enr = nw_annot.enrich(gene_set='start_genes')
NetAn.plot_utils.plot_enrichments(MS_enr, sorton='pval', plot_Ntop=10)

MS_enr.to_csv('MS_enrichments.csv')
MS_cluster_exp['cluster_0'][0].to_csv('cluster_1.csv')
MS_cluster_exp['cluster_1'][0].to_csv('cluster_2.csv')


nw_annot.add_close_genes(augment_factor=1)
nw_annot.enrich(gene_set='genes_exp')
MS_enr_exp = nw_annot.enrich(gene_set='genes_exp')

MS_enr_exp.to_csv('MS_enrichments_expanded.csv')

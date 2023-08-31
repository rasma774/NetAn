import pandas as pd
import NetAn

genes = pd.read_csv('netGO_genes.csv').values.T[0]

df = pd.read_csv('/home/rasmus/Documents/projects/nw_annot/data/netGO_benchmark/netGO_clean.csv',
                 index_col=0,
                 header=None)



nt = NetAn.NetAn(genes)
nt.add_close_genes()

df_connected = nt.enrich(annotations=df, gene_set='all', connected_components=True)
df_expanded = nt.enrich(annotations=df, gene_set='all')
df_cliques = nt.enrich(annotations=df, gene_set='all', use_cliques=True)
df_cluster = nt.enrich(annotations=df, gene_set='all', n_clusters=2)


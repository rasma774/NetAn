import pandas as pd
import numpy as np
import sklearn.cluster
import networkx as nx
import scipy.stats as sts
import os

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2023 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

PATH = os.path.dirname(os.path.abspath(__file__)).strip('utils')
np.random.seed(0)

def load_nw(nw, cutoff=0):
    if nw == 'STRINGdb':
        df = pd.read_csv(PATH + 'data/STRINGdb/clean_string.csv')
        df = df.set_index('p1')
        df = df[df.score >= cutoff]
    return df

def map2nw(start_nodes, nw, augment_factor=2, norm=True):
    if type(start_nodes) is list:
        start_nodes = np.array(start_nodes)
    elif type(start_nodes) is not np.ndarray:
        raise ValueError('start_nodes must be np.array or list')

    in_nw = np.in1d(start_nodes, np.unique(nw.index))
    if not np.all(in_nw):
        print(str(np.sum(~in_nw)) +
              ' of ' +
              str(len(in_nw)) +
              ' genes are not in network')

    start_nodes_filtered = start_nodes[in_nw]
    score_sum = nw.loc[start_nodes_filtered].groupby(nw.columns[0]).sum()

    if norm:
        norm_vals = pd.read_csv(PATH +
                                'data/nw_backgrounds/STRINGdb_background.csv',
                                index_col=0)

        score_sum = score_sum/norm_vals.reindex(score_sum.index)

    score_sum = score_sum.sort_values('score', ascending=False)
    score_sum = score_sum[~score_sum.index.isin(start_nodes)]
    gene_exp = score_sum.index[:int(augment_factor*len(start_nodes))]
    return np.array(gene_exp)


def cluster_nw(genes, nw, n_clusters=None):
    genes = genes[np.in1d(genes, np.unique(nw.index))]
    subnw = nw.loc[genes]
    subnw = subnw[subnw.iloc[:, 0].isin(genes)]
    subnw = subnw.pivot(columns=nw.columns[0])
    subnw = subnw.fillna(0)

    if n_clusters is None:
        kmeans = sklearn.cluster.KMeans(n_init='auto')
    else:
        kmeans = sklearn.cluster.KMeans(n_clusters=n_clusters)

    kmeans.fit(subnw.values)
    return kmeans.labels_, subnw.index


def filter_walk(walked_genes, cutoff=2, test=None):
    # possibly add other tests in future releases
    gene_name, count = np.unique(walked_genes, return_counts=True)
    gene_name = gene_name[count > cutoff]
    return gene_name


def _build_graph(nw):
    print('Building graph in networkx')
    graph = nx.Graph(list(zip(nw.index, nw.iloc[:, 0])))
    nx.connected_components(graph)

    # We extract the biggest component
    components = [graph.subgraph(comp).copy() for comp in nx.connected_components(graph)]
    lens = [len(c.nodes) for c in components]
    biggest_component = components[np.argmax(lens)]

    print('Lost ', 1 - (len(biggest_component.nodes)/len(graph.nodes)),
          'of genes by extracting largest connected component of the network')
    return biggest_component

def _get_dists(network, gene_list1, gene_list2, n_rand=100):
    # This is an internal function that runs a number of random distance calcs
    rand_lens = []

    randpos_1 = np.random.randint(0, len(gene_list1), (n_rand, 1))
    randpos_2 = np.random.randint(0, len(gene_list2), (n_rand, 1))

    genes1 = np.array(gene_list1)[randpos_1[:, 0]]
    genes2 = np.array(gene_list2)[randpos_2[:, 0]]

    for i in range(len(randpos_1)):
        if genes1[i] == genes2[i]:
            continue

        rand_lens.append(nx.shortest_path_length(network,
                                                  source=genes1[i],
                                                  target=genes2[i],
                                                  )
                          )
    return np.array(rand_lens)

def _cluster_dist_until_convergence(network, genes1, genes2):
    # This internal function calculates the average shortest pathway by random
    # sampling until convergence
    lens = _get_dists(network, genes1, genes2)
    q = np.inf
    while not (0.99 < q < 1.01):
        print(q)
        lens_new = _get_dists(network, genes1, genes1)
        new_lens = np.concatenate((lens, lens_new))
        q = lens.mean()/new_lens.mean()
        lens = new_lens
    return new_lens


def compare_distance(cluster_res, nw, max_size=300):
    # if size of cluster is above 300, then monte-carlo estimate until converge
    network = _build_graph(nw.copy())

    lens = _get_dists(network, network.nodes, network.nodes)

    q = np.inf
    while not (0.99 < q < 1.01):
        lens_new = _get_dists(network, network.nodes, network.nodes)
        new_lens = np.concatenate((lens, lens_new))
        q = lens.mean()/new_lens.mean()
        lens = new_lens
    mean_background = new_lens.mean()

    keys = list(cluster_res.keys())
    res_mu = np.zeros((len(keys), len(keys)))
    res_mu[:] = np.nan
    p_mannwhitney = res_mu.copy()
    for i in range(len(keys)):
        for j in range(i, len(keys)):
            print(i,j)
            genes1 = cluster_res[keys[i]][1]
            genes2 = cluster_res[keys[j]][1]

            if max(len(genes1), len(genes2)) > max_size:
                lens_tmp = _cluster_dist_until_convergence(network,
                                                           genes1,
                                                           genes2)
            else:
                lens_tmp = []
                for g1 in genes1:
                    for g2 in genes2:
                        lens_tmp.append(
                            nx.shortest_path_length(network,
                                                    source=g1,
                                                    target=g2,
                                                    ))
            res_mu[i, j] = np.mean(lens_tmp)
            res_mu[j, i] = res_mu[i, j]
            p_mannwhitney[i, j] =  sts.mannwhitneyu(new_lens, lens_tmp)[1]
            p_mannwhitney[j, i] = p_mannwhitney[i, j]
    mean_fold = res_mu/mean_background
    return mean_fold, p_mannwhitney

def get_connected_components(genes, nw):
    nw = nw.copy()
    nw = nw[nw.index.isin(genes)]
    nw = nw[nw.iloc[:, 0].isin(genes)]

    print('Building graph in networkx')
    graph = nx.Graph(list(zip(nw.index, nw.iloc[:, 0])))
    nx.connected_components(graph)
    components = [graph.subgraph(comp).copy() for comp in nx.connected_components(graph)]

    return [np.array(c.nodes) for c in components]

def get_cliques(genes, nw):
    nw = nw.copy()
    nw = nw[nw.index.isin(genes)]
    nw = nw[nw.iloc[:, 0].isin(genes)]

    print('Building graph in networkx')
    graph = nx.Graph(list(zip(nw.index, nw.iloc[:, 0])))

    return [np.array(clq) for clq in nx.find_cliques(graph)]

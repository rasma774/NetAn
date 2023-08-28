import utils.nw_utils as nw_utils
import utils.enrich_utils as enrich_utils
import utils.plot_utils as plot_utils
import numpy as np

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2023 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

class NetAn:
    def __init__(self, genes, network='STRINGdb', nw_cutoff=0):
        """
        Parameters
        ----------
        genes : list or numpy array
            Should contain all genes to be tested
        network : str or pandas DataFrame, optional
            The default is 'STRINGdb', which maps to the built-in
            protein-protein interaction network STRING. Alternatively, you can
            include a custom network as a pandas DataFrame. More information
            can be found on wwww.github.com/XXX
        nw_cutoff : float, optional
            A cutoff for the STRING network confidence. Between 0 and 1.
            The default is 0.

        Notes
        -----

        This is the Network Annotation Enrichment package, NetAn, which takes a
        list of genes and uses network-based approaches such as network
        clustering and inference of closely related genes to include local
        neighbours.

        Example Usage
        -------------
        >>> from NetAn import NetAn
        >>> nt = NetAn(['GATA3', 'STAT4', 'IFG3']) # Example using three genes
        >>> nt.add_close_genes() # Expand the gene set for later enrichment
        >>> _ = nt.genes_exp # The added genes are stored here
        >>> # Test for annotation enr in input and expanded genes
        >>> enrichments = nt.enrich(gene_set='all',
                                    n_clusters=2,
                                    )
        >>> f, ax = plot_enrichment_clusters()


        """
        if network == 'STRINGdb':
            self.network = nw_utils.load_nw('STRINGdb', cutoff=nw_cutoff)
        self.genes = np.array(genes)
        self.genes_exp = None
        self.genes_all = None
        self.clusters = None


    def add_close_genes(self, augment_factor=2, default_norm=True):
        """
        Identify genes that are close to the input genes, based on the input
        network. The new genes are saved as the attribute NetAn.genes_exp
        and the total gene set is stored as NetAn.genes_all

        Parameters
        ----------
        augment_factor : int, optional
            The fold of number of genes to be inferred, as compared to number of
            genes in the input. Example: an augment_factor of 2 and 10 input
            genes would infer 20 new genes. The default value of augment_factor
            is 2.
        norm : bool, optional
            Select whether of not the network should be normed. Only implemented
            when NetAn uses the default STRINGdb network. The default is True.

        Returns
        -------
        None.

        Example usage
        -------------
        >>> # Assuming an NetAn obj called nt
        >>> nt.add_close_genes(augment_factor=1)
        """
        self.genes_exp = nw_utils.map2nw(self.genes,
                                         self.network,
                                         augment_factor=augment_factor,
                                         norm=default_norm)
        self.genes_all = np.unique(
            np.concatenate([self.genes_exp, self.genes])
                           )


    def enrich(self,
               annotations='GO',
               gene_set='start_genes',
               n_clusters=-1,
               connected_components=False,
               use_cliques=False,
               clique_size_thresh=5,
               nw_thresh=0,
               ):
        """
        Identify significant annotation enrichments in either gene network
        clusters, or in the gene sets.

        Parameters
        ----------
        annotations : string or pandas DataFrame
            Set what annotation-gene sets to include in the study. The default
            is the accompanying GO set.
        gene_set : str. 'start_genes' | 'genes_exp' | 'all', optional
            Set what genes to use in the over-expression analysis. If the
            NetAn.add_close_genes method has been run, the 'genes_exp' can be
            used to analyse the added neighbours of the start genes, or 'all'
            can be used to sets both the start genes and the added ones.
            The default is 'start_genes'.
        n_cluster : int, optional
            Select whether to cluster the genes. Clustering is based on the
            network, with n_clusters setting the number of used clusters the
            genes should be divided into. The default is -1, meaning no
            clustering.
        connected_components : bool
            When set to True, NetAn overrides the n_cluster parameter and the
            cliques parameter and extracts the connected components when the
            genes are mapped to the network, and does an overrepresentation
            analysis for all components. The default is False.
        use_cliques : bool
            Extracts the maximal cliques from the input genes when mapped to
            the network, and does an overrepresentation analysis for each
            clique. The default is False.
        clique_size_thresh : int
            Sets a threshold for how big the cliques should be to be the
            included in the overrepresentation analysis. The default is 5.
        nw_thresh : float; 0 < nw_thresh < 1 for built-in network
            Set a threshold for the network when used to subdivide the genes
            into subsets (i.e. into connected components, cliques, clustering)
            The default is 0, meaning no filtering

        Raises
        ------
        ValueError

        Returns
        -------
        cluster_res : dict or pandas DataFrame
            If clustering is used, cluster_res is a dict of the results of each
            cluster, and the genes in it. Otherwise, it is a DataFrame with the
            results.

        """
        # Some initial set ups
        if gene_set == 'genes_exp':
            genes = self.genes_exp
        elif gene_set == 'start_genes':
            genes = self.genes
        elif gene_set == 'all':
            genes = self.genes_all
        else:
            genes = gene_set

        if genes is None:
            raise ValueError('Genes not set')

        nw_tmp = self.network[self.network.iloc[:, -1] > nw_thresh]

        if connected_components:
            con_comp = nw_utils.get_connected_components(genes, nw_tmp)

            # This is if there is only one connected component
            if not type(con_comp[0]) is np.ndarray:
                con_comp = [con_comp]
            comp_res = {}
            for i, compgenes in enumerate(con_comp):
                comp_res['component_' + str(i)] = [enrich_utils.calc_enrich(
                    compgenes,
                    annotations),
                    compgenes]
            return comp_res

        elif use_cliques:

            cliques = nw_utils.get_cliques(genes, nw_tmp)
            comp_res = {}
            for i, compgenes in enumerate(cliques):
                if len(compgenes) < clique_size_thresh:
                    continue
                comp_res['component_' + str(i)] = [enrich_utils.calc_enrich(
                    compgenes,
                    annotations),
                    compgenes]
            return comp_res


        elif n_clusters > 0:
            clusters, cluster_genes = nw_utils.cluster_nw(genes,
                                                          nw_tmp,
                                                          n_clusters=n_clusters)

            cluster_res = {}
            for x in np.unique(clusters):
                tmp_genes = cluster_genes[clusters == x]
                cluster_res['cluster_' + str(x)] = [
                    enrich_utils.calc_enrich(tmp_genes, annotations),
                    tmp_genes,
                    ]

            self.clusters = cluster_res

            return cluster_res

        return enrich_utils.calc_enrich(genes, annotations)

    def get_cluster_distances(self):
        """
        Can be run after the NetAn.enrich method has been applied. Calculates
        a distance matrix between the clusters

        Returns
        -------
        dists : numpy.ndarray
            contains the normed distances between the clusters, compared to the
            distances of the full network.
        pvals : numpy.ndarray
            The p-values for H0: the distances have the same mean as between two
            random paths in the nework.

        Notes
        -----
        >>> # Run example, assuming that the NetAn.enrich with n_clusters > 1
        >>> (here 2) as an input param has been run.
        >>> dists, pvals = nt.get_cluster_distances()

        >>> print(dists)
                [[0.53, 0.9],
                 [0.9, 0.23]]

            In this example, there is a quite long distance between the two
            clusters (0.9 of average), but the internal distances are small,
            meaning that the two clusters were well defined. Cluster 2 is
            especially close and well defined, as seen by the low number 0.23

        """
        if self.clusters is None:
            raise ValueError('No clusters created')

        dists, pvals = nw_utils.compare_distance(self.clusters, self.network)
        return dists, pvals

    def cite(self):
        print('Please support this software by citing XXX.')
        print('Moreover, if the STRINGdb material has been used, also cite')
        print('https://string-db.org/cgi/about?footer_active_subpage=references')

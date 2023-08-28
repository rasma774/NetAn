import pandas as pd
import scipy.stats as sts
import numpy as np
from statsmodels.stats.multitest import multipletests
import os

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2023 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'


PATH = os.path.dirname(os.path.abspath(__file__)).strip('utils')


def calc_enrich(genes, annotations='GO', ngenes_background=22000, thresh=20):
    if type(annotations) is str:
        if annotations == 'GO':
            annotations = pd.read_csv(PATH + 'data/GO/GO_clean.csv',
                                      index_col=0)


    unique_annots, counts = np.unique(annotations.index, return_counts=True)
    unique_annots = unique_annots[counts > thresh]
    annotations = annotations.loc[unique_annots]

    all_res = {}
    print('Calculating Fisher test')
    for term in np.unique(annotations.index):
        # Do a fisher test
        annot_genes = annotations.loc[term].values.T[0]
        annot_genes = np.unique(annot_genes)

        # Now we prepare for the fisher test
        in_both = np.in1d(genes, annot_genes).sum()
        only_genes = len(genes) - in_both
        only_annot = len(annot_genes) - in_both
        not_in_any = ngenes_background - (in_both + only_annot + only_genes)

        fishres = sts.fisher_exact([[in_both, only_annot],
                                    [only_genes, not_in_any]],
                                   alternative='greater')
        all_res[term] = {'OR': fishres[0], 'pval':fishres[1]}

    df = pd.DataFrame(all_res).transpose()

    fdr = multipletests(df.pval.values, method='fdr_bh')
    df['p_adj'] = fdr[1]

    return df.sort_values('pval')

import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt

incs = []
ORs = []
Ps = []
all_ORs = []
all_new_ORs_in_old = []
for n_aug in [.5, 1, 1.5, 2, 2.5]:
    with open('power_increase_nAug_' + str(n_aug) + '.pkl', 'rb') as f:
        res = pkl.load(f)
    n_tmp = np.array(res['n_in_sets'])
    n_tmp = np.array(res['n_in_sets'])
    fold_tmp = n_tmp[:, 2]/n_tmp[:, 1]
    fold_tmp = fold_tmp[~np.isinf(fold_tmp)]
    incs.append(np.nanmean(fold_tmp))

    Ps.append(np.nanmedian(np.array(res['only2_p'])[:, 0]))
    ORs.append(np.nanmedian(np.array(res['only2_OR'])[:, 0]))
    all_ORs.append(np.array(res['only2_OR'])[:, 0])
    all_new_ORs_in_old.append(np.array(res['only1_OR'])[:, 1])

for i in range(len(all_ORs)):
    all_ORs[i] = all_ORs[i][~np.isnan(all_ORs[i])]
    all_new_ORs_in_old[i] = all_new_ORs_in_old[i][~np.isnan(all_new_ORs_in_old[i])]

f, ax = plt.subplots(1, figsize=(4 , 4))
ax.boxplot(all_ORs, showfliers=False, boxprops={'linewidth':1.3}, whiskerprops={'linewidth':1.3})#ax.plot(np.nanmedian(all_ORs, 1), 'x', markersize=7)
ax.set_xlim([0,6.5])
f.savefig('../res/ORs_not_sign_in_1st_line_boxplot.svg')


f, ax = plt.subplots(1, figsize=(4 , 4))
#_corr_ax(ax)
ax.boxplot(all_new_ORs_in_old, showfliers=False, boxprops={'linewidth':1.3}, whiskerprops={'linewidth':1.3})#ax.plot(np.nanmedian(all_ORs, 1), 'x', markersize=7)
ax.set_xlim([0,6.5])
f.savefig('../res/ORs_not_sign_in_2nd_line_boxplot.svg')

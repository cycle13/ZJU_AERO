'''
Description: test gamma 
Author: Hejun Xie
Date: 2020-12-15 09:11:39
LastEditors: Hejun Xie
LastEditTime: 2020-12-30 23:03:43
'''

A_AR_LAMBDA_AGG = 8.42
B_AR_LAMBDA_AGG = -0.57
A_AR_M_AGG =  0.053
B_AR_M_AGG = 0.79

# A_AR_LAMBDA_AGG = 4.0
# B_AR_LAMBDA_AGG = -0.50
# A_AR_M_AGG =  0.04
# B_AR_M_AGG =  0.70

from scipy import stats
import numpy as np

if __name__ == "__main__":
    D_list = np.arange(2, 15., 0.5)
    asp = np.arange(0.5, 3.5, 0.01)

    X, Y = np.meshgrid(D_list, asp)

    ar_lambda = A_AR_LAMBDA_AGG * D_list ** B_AR_LAMBDA_AGG
    ar_mu = A_AR_M_AGG * D_list ** B_AR_M_AGG
    ar_loc = np.ones((len(D_list),))

    result = np.zeros((len(D_list), len(asp)), dtype='float32')

    for l in zip(range(len(D_list)), ar_lambda, ar_loc, ar_mu):
        gamm = stats.gamma(l[1], l[2], l[3])
        wei = gamm.pdf(asp)
        # wei /= np.sum(wei) # renormalization
        result[l[0], :] = wei

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'

    fig, ax = plt.subplots(figsize=(10,8))
    cs = ax.contourf(X, Y, result.T, cmap='RdBu_r', levels=30)
    cbar = fig.colorbar(cs)

    cbar.ax.set_ylabel('Probability Distribution', fontsize=14)

    ax.set_xlabel('Dmax [mm]', fontsize=14)
    ax.set_ylabel('Aspect Ratio', fontsize=14)

    plt.savefig('test_gamma.png', dpi=300)

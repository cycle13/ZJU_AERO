'''
Description: test particle size distribution for 
different microhphysics schemes
Author: Hejun Xie
Date: 2021-10-18 15:57:32
LastEditors: Hejun Xie
LastEditTime: 2021-10-18 20:45:33
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np

# Local imports
from ZJU_AERO.config.cfg import createConfig
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.hydro.hydrometeor import Rain, Graupel

if __name__ == "__main__":
    createConfig('./option_files/thompson_test.yml')
    constants.update()
    
    r_wsm6 = Rain('1mom', scheme_name='wsm6')
    g_wsm6 = Graupel('1mom', scheme_name='wsm6')

    r_thp = Rain('1mom', scheme_name='thompson')
    g_thp = Graupel('1mom', scheme_name='thompson')

    Dr = r_wsm6.list_D
    Dg = g_wsm6.list_D

    lwcs = [1E-5, 1E-4, 1E-3, 1E-2]  # [kg m-3]

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'
    mpl.rcParams['text.usetex'] = True

    lwcs = [1E-3, 2E-3, 5E-3, 1E-2]  # [kg m-3]
    text = [r"LWC = 1 [$ g \cdot m^{-3} $]", r"LWC = 2 [$ g \cdot m^{-3} $]",
            r"LWC = 5 [$ g \cdot m^{-3} $]", r"LWC = 10 [$ g \cdot m^{-3} $]"]

    # # size distribution
    fig, axes = plt.subplots(2, 2, figsize=(10, 15), sharey=True, sharex=True)

    for ilwc in range(len(lwcs)):
        ax = axes[ilwc // 2, ilwc % 2]

        r_wsm6.set_psd(lwcs[ilwc])
        nd_wsm6 = r_wsm6.get_N(Dr).squeeze() # [mm-1 m-3]
        mass_wsm6 = r_wsm6.a * Dr**r_wsm6.b * nd_wsm6 # [kg mm-1 m-3]

        r_thp.set_psd(lwcs[ilwc])
        nd_thp = r_thp.get_N(Dr).squeeze() # [mm-1 m-3]
        mass_thp = r_thp.a * Dr**r_thp.b * nd_thp # [kg mm-1 m-3]
        
        ax.hist(Dr, len(Dr), weights=mass_wsm6 * 1E3, histtype='stepfilled', density=False,
            label='wsm6', facecolor='r', alpha=0.5)
        
        ax.hist(Dr, len(Dr), weights=mass_thp * 1E3, histtype='stepfilled', density=False,
            label='thompson', facecolor='b', alpha=0.5)

        ax.set_xscale('log')

        ax.text(0.25, 2, text[ilwc], size=14, ha="center", va="center")

        ax.set_ylabel(r"Mass Distribution M(D) [$ g \cdot m^{-3} \cdot mm^{-1} $]", fontsize=14)
        ax.set_xlabel(r"Maximum Dimension $D_{max}$ [$ mm $]", fontsize=16)

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

    
    axes[0, 0].legend(loc='best', fontsize=15, frameon=False)
    fig.tight_layout(rect=[0.08, 0.35, 0.97, 0.97])
    plt.subplots_adjust(left=0.08, bottom=0.35, right=0.97, top=0.97)
    aux_ax = plt.axes([0.07, 0.07, 0.9, 0.2])

    QM = np.logspace(-5, -2, 1000)

    r_wsm6.set_N0(QM)
    r_thp.set_N0(QM)
    N0_wsm6 = r_wsm6.N0 + np.zeros(QM.shape, dtype='float32')
    N0_thp  = r_thp.N0 + np.zeros(QM.shape, dtype='float32')
    
    aux_ax.plot(QM*1E3, N0_wsm6, label='wsm6', color='r')
    aux_ax.plot(QM*1E3, N0_thp, label='thompson', color='b')

    aux_ax.plot(QM[QM<r_thp.Q0]*1E3, r_thp.N1 + np.zeros(QM[QM<r_thp.Q0].shape, dtype='float32'), 'k--')
    aux_ax.plot(QM[QM>r_thp.Q0]*1E3, r_thp.N2 + np.zeros(QM[QM>r_thp.Q0].shape, dtype='float32'), 'k--')
    plt.axvline(x=r_thp.Q0*1E3, color='k', ls='-')
    plt.axhline(y=N0_thp[np.argmin(np.abs(QM-r_thp.Q0))], color='k', ls='-')

    eq1 = r"$N_{0}=\frac{N_1-N_2}{2} \cdot \tanh(\frac{q_{r0}-q_r}{4q_{r0}}) + \frac{N_1+N_2}{2}$"
    aux_ax.text(3E-2, 2E5, eq1, color="k", fontsize=16,
        horizontalalignment="center", verticalalignment="center")
    
    aux_ax.set_xscale('log')
    aux_ax.set_yscale('log')

    aux_ax.set_ylabel(r"Intercept parameter $N_{0}$ [$ mm^{-1} \cdot m^{-3} $]", fontsize=14)
    aux_ax.set_xlabel(r"Liquid Water Content [$ g \cdot m^{-3} $]", fontsize=14)

    aux_ax.legend(loc='best', fontsize=14, frameon=False)

    aux_ax.spines['bottom'].set_linewidth(1.5)
    aux_ax.spines['left'].set_linewidth(1.5)
    aux_ax.spines['right'].set_linewidth(1.5)
    aux_ax.spines['top'].set_linewidth(1.5)
    
    plt.savefig('test_psd_rain.png', dpi=300)
    plt.close()


    # # size distribution
    fig, axes = plt.subplots(2, 2, figsize=(10, 15), sharey=True, sharex=True)

    iwcs = [5E-4, 1E-3, 2E-3, 5E-3]  # [kg m-3]
    text = [r"IWC = 0.5 [$ g \cdot m^{-3} $]", r"IWC = 1 [$ g \cdot m^{-3} $]",
            r"IWC = 2 [$ g \cdot m^{-3} $]", r"IWC = 5 [$ g \cdot m^{-3} $]"]

    for iiwc in range(len(iwcs)):
        ax = axes[iiwc // 2, iiwc % 2]

        g_wsm6.set_psd(iwcs[iiwc])
        nd_wsm6 = g_wsm6.get_N(Dg).squeeze() # [mm-1 m-3]
        mass_wsm6 = g_wsm6.a * Dg**g_wsm6.b * nd_wsm6 # [kg mm-1 m-3]

        g_thp.set_psd(iwcs[iiwc])
        nd_thp = g_thp.get_N(Dg).squeeze() # [mm-1 m-3]
        mass_thp = g_thp.a * Dg**g_thp.b * nd_thp # [kg mm-1 m-3]
        
        ax.hist(Dg, len(Dg), weights=mass_wsm6 * 1E3, histtype='stepfilled', density=False,
            label='wsm6', facecolor='r', alpha=0.5)
        
        ax.hist(Dg, len(Dg), weights=mass_thp * 1E3, histtype='stepfilled', density=False,
            label='thompson', facecolor='b', alpha=0.5)

        ax.set_xscale('log')

        ax.text(0.5, 0.5, text[iiwc], size=14, ha="center", va="center")

        ax.set_ylabel(r"Mass Distribution M(D) [$ g \cdot m^{-3} \cdot mm^{-1} $]", fontsize=14)
        ax.set_xlabel(r"Maximum Dimension $D_{max}$ [$ mm $]", fontsize=14)

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

    
    axes[0, 0].legend(loc='best', fontsize=14, frameon=False)
    
    fig.tight_layout(rect=[0.08, 0.35, 0.97, 0.97])
    plt.subplots_adjust(left=0.08, bottom=0.35, right=0.97, top=0.97)
    aux_ax = plt.axes([0.07, 0.07, 0.9, 0.2])

    QM = np.logspace(-5, -2, 1000)

    g_wsm6.set_N0(QM)
    g_thp.set_N0(QM)
    N0_wsm6 = g_wsm6.N0 + np.zeros(QM.shape, dtype='float32')
    N0_thp  = g_thp.N0 + np.zeros(QM.shape, dtype='float32')
    
    aux_ax.plot(QM*1E3, N0_wsm6, label='wsm6', color='r')
    aux_ax.plot(QM*1E3, N0_thp, label='thompson', color='b')

    aux_ax.set_xscale('log')
    aux_ax.set_yscale('log')

    aux_ax.set_ylabel(r"Intercept parameter $N_{0}$ [$ mm^{-1} \cdot m^{-3} $]", fontsize=14)
    aux_ax.set_xlabel(r"Ice Water Content [$ g \cdot m^{-3} $]", fontsize=14)

    eq2 = r"$N_{0}=max\left(N_1,  min\left( \frac{\xi}{q_g}, N_2 \right) \right)$"
    aux_ax.text(3E-2, 2E2, eq2, color="k", fontsize=18,
        horizontalalignment="center", verticalalignment="center")

    aux_ax.legend(loc='best', fontsize=14, frameon=False)

    aux_ax.spines['bottom'].set_linewidth(1.5)
    aux_ax.spines['left'].set_linewidth(1.5)
    aux_ax.spines['right'].set_linewidth(1.5)
    aux_ax.spines['top'].set_linewidth(1.5)

    
    plt.savefig('test_psd_graupel.png', dpi=300)
    plt.close()


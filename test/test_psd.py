'''
Description: test particle size distribution for 
different microhphysics schemes
Author: Hejun Xie
Date: 2021-10-18 15:57:32
LastEditors: Hejun Xie
LastEditTime: 2021-10-18 17:50:55
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

    lwcs = [1E-3, 2E-3, 5E-3, 1E-2]  # [kg m-3]
    text = [r"LWC = 1 [$ g \cdot m^{-3} $]", r"LWC = 2 [$ g \cdot m^{-3} $]",
            r"LWC = 5 [$ g \cdot m^{-3} $]", r"LWC = 10 [$ g \cdot m^{-3} $]"]

    # # size distribution
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey=True, sharex=True)

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

        ax.text(0.25, 2, text[ilwc], size=12, ha="center", va="center")

        ax.set_ylabel(r"Mass Distribution M(D) [$ g \cdot m^{-3} \cdot mm^{-1} $]", fontsize=12)
        ax.set_xlabel(r"Maximum Dimension $D_{max}$ [$ mm $]", fontsize=12)

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

    
    axes[0, 0].legend(loc='best', fontsize=12, frameon=False)
    fig.tight_layout()
    plt.savefig('test_psd_rain.png', dpi=300)
    plt.close()


    # # size distribution
    fig, axes = plt.subplots(2, 2, figsize=(10, 10), sharey=True, sharex=True)

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

        ax.text(0.5, 0.5, text[iiwc], size=12, ha="center", va="center")

        ax.set_ylabel(r"Mass Distribution M(D) [$ g \cdot m^{-3} \cdot mm^{-1} $]", fontsize=12)
        ax.set_xlabel(r"Maximum Dimension $D_{max}$ [$ mm $]", fontsize=12)

        ax.spines['bottom'].set_linewidth(1.5)
        ax.spines['left'].set_linewidth(1.5)
        ax.spines['right'].set_linewidth(1.5)
        ax.spines['top'].set_linewidth(1.5)

    
    axes[0, 0].legend(loc='best', fontsize=12, frameon=False)
    fig.tight_layout()
    plt.savefig('test_psd_graupel.png', dpi=300)
    plt.close()


'''
Description: test terminal velocity for different microphysics schemes
Author: Hejun Xie
Date: 2021-10-18 14:57:39
LastEditors: Hejun Xie
LastEditTime: 2021-10-18 15:17:41
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

# Global imports
import numpy as np

# Local imports
from ZJU_AERO.config.cfg import createConfig
from ZJU_AERO.const import global_constants as constants
from ZJU_AERO.hydro.hydrometeor import Snow, Rain, Graupel

if __name__ == "__main__":

    createConfig('./option_files/thompson_test.yml')
    constants.update()

    Ds = np.linspace(0.2, 20, 1024)
    Dr = np.linspace(0.1, 15, 1024)
    Dg = np.linspace(0.2, 15, 1024)

    s = Snow('1mom', scheme_name='wsm6')
    r = Rain('1mom', scheme_name='wsm6')
    g = Graupel('1mom', scheme_name='wsm6')

    s_thp = Snow('1mom', scheme_name='thompson')
    r_thp = Rain('1mom', scheme_name='thompson')
    g_thp = Graupel('1mom', scheme_name='thompson')

    vt_s = s.get_V(Ds)
    vt_r = r.get_V(Dr)
    vt_g = g.get_V(Dg)

    vt_s_thp = s_thp.get_V(Ds)
    vt_r_thp = r_thp.get_V(Dr)
    vt_g_thp = g_thp.get_V(Dg)

    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'
    
    fig, axes = plt.subplots(3, 1, figsize=(10,12), sharex=False)
    
    axes[0].plot(Ds, vt_s, color='k', label='wsm6', ls='-')
    axes[0].plot(Ds, vt_s_thp, color='k', label='thompson', ls='--')

    axes[0].set_ylabel(r'Snow Terminal velocity $v_{t}$ [m/s]', fontsize=14)
    axes[0].set_xlabel(r'Snow Maximum Diemnsion $D_{max}$ [m/s]', fontsize=14)

    axes[1].plot(Dr, vt_r, color='k', label='wsm6', ls='-')
    axes[1].plot(Dr, vt_r_thp, color='k', label='thompson', ls='--')

    axes[1].set_ylabel(r'Rain Terminal velocity $v_{t}$ [m/s]', fontsize=14)
    axes[1].set_xlabel(r'Rain Maximum Diemnsion $D_{max}$ [m/s]', fontsize=14)

    axes[2].plot(Dg, vt_g, color='k', label='wsm6', ls='-')
    axes[2].plot(Dg, vt_g_thp, color='k', label='thompson', ls='--')

    axes[2].set_ylabel(r'Graupel Terminal velocity $v_{t}$ [m/s]', fontsize=14)
    axes[2].set_xlabel(r'Graupel Maximum Diemnsion $D_{max}$ [m/s]', fontsize=14)

    axes[0].legend(loc='best', frameon=False, fontsize=12)
    fig.tight_layout()
    plt.savefig('test_vt.png', dpi=300)
    plt.close()
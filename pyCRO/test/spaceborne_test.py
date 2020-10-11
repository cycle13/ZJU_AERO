'''
Description: test for spaceborne radar
Author: Hejun Xie
Date: 2020-10-10 10:44:15
LastEditors: Hejun Xie
LastEditTime: 2020-10-11 17:15:20
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/pyCRO/')

# Global imports
import numpy as np
import pickle
# import pyart

# Local imports
import pyCRO
import pyart

LOAD_MODEL = True
LOAD_RADAR = True
DEG = r'$^\circ$'

cmap = {'ZH':'pyart_Carbone11', 'RVEL': 'pyart_BuOr8', 'ZDR': 'pyart_Carbone17',
'KDP': 'pyart_EWilson17', 'PHIDP': 'pyart_Carbone42', 'RHOHV': 'pyart_GrMg16'}
clevels = {'ZH':np.arange(0,55,5), 'ZDR':np.linspace(0., 0.1, 50)}
units = {'ZH':'[dBZ]', 'ZDR':'[dBZ]'}

if __name__ == "__main__":

    if not LOAD_RADAR:
        FILENAME = '../../../cosmo_pol/pathos/WRF/wsm6_test/wrfout_d01_2019-05-17_00_00_00'
        a = pyCRO.RadarOperator(options_file='./option_files/spaceborne_test.yml')
        a.load_model_file(FILENAME, itime=40, load_pickle=LOAD_MODEL, pickle_file='mdl.pkl')

    # print(a.dic_vars['T'])

    if not LOAD_RADAR:
        r = a.get_spaceborne_swath('test')
        with open("./swath.pkl", "wb") as f:
            pickle.dump(r, f)
    else:
        with open("./swath.pkl", "rb") as f:
            r = pickle.load(f)
    
    # exit()
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = 'serif'
    from mpl_toolkits.basemap import Basemap

    var = 'ZH'

    fig, ax = plt.subplots(figsize=(10,8))

    slat = 38.0
    elat = 42.0
    slon = 115.0
    elon = 119.0

    map = Basemap(projection='cyl', llcrnrlat=slat, urcrnrlat=elat, llcrnrlon=slon, urcrnrlon=elon, resolution='h', ax=ax)
    map.drawparallels(np.arange(slat, elat+1, 1.0), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
    map.drawmeridians(np.arange(slon, elon+1, 1.0), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])
    map.drawcoastlines()

    x = r.lons
    y = r.lats
    h = r.heights
    z = 10*np.log10(r.data[var])

    zlevel = 35
    
    im = map.pcolormesh(x[...,zlevel], y[...,zlevel], z[...,zlevel], cmap=cmap[var])

    ax.set_title('Spaceborne Radar ' + var, fontsize=16)
    
    cb = fig.colorbar(im, ax=ax)
    cb.ax.set_ylabel('Horizontal Reflectivity [dBZ]', fontsize=14)

    plt.savefig('swath.png',dpi=300, bbox_inches='tight')
    exit()

    # start swath 3D plot
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import axes3d
    
    var = 'ZH'

    fig = plt.figure()
    ax = fig.gca(projection='3d')

    x = r.lons
    y = r.lats
    h = r.heights
    z = 10*np.log10(r.data[var])
    # levels = clevels[var]
    zgrids = [5, 15, 25, 35]

    csets = []
    for zgrid in zgrids:
        levels = np.arange(np.nanmin(z[:,:,zgrid]), np.nanmax(z[:,:,zgrid]), 1)
        csets.append(ax.contourf(x[:,:,zgrid], y[:,:,zgrid], z[:,:,zgrid],
            levels, offset=h[0,0,zgrid]/1000., cmap = cmap[var]))
        
    cbar = fig.colorbar(csets[0])
    cbar.ax.set_ylabel(units[var])

    ax.set_xlim(115, 117.5)
    ax.set_ylim(38, 41)
    ax.set_zlim(0, 10)

    ax.set_xlabel('longitude')
    ax.set_ylabel('latitude')
    ax.set_zlabel('heights [km]')
    ax.set_title('Spaceborne Radar ' + var)

    plt.show()

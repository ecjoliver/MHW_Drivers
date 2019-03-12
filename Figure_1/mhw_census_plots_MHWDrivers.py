'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to 
  select SST time series around the globe

'''

import numpy as np
import scipy.signal as sig
from scipy import linalg
from scipy import stats

import ecoliver as ecj

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm
#import cartopy.crs as ccrs


#
# Load data and make plots
#

outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016'
#outfile = '/media/ecoliver/DataOne/mokami_move_backup/data/MHWs/Trends/mhw_census.2016'

data = np.load(outfile+'.npz')
lon_map = data['lon_map']
lat_map = data['lat_map']
SST_mean = data['SST_mean']
MHW_total = data['MHW_total']
MHW_cnt = data['MHW_cnt']
MHW_dur = data['MHW_dur']
MHW_max = data['MHW_max']
MHW_mean = data['MHW_mean']
MHW_cum = data['MHW_cum']
MHW_td = data['MHW_td']
MHW_tc = data['MHW_tc']

SST_mean[SST_mean==0] = np.nan
MHW_cnt[MHW_cnt==0] = np.nan
MHW_dur[MHW_cnt==0] = np.nan
MHW_max[MHW_cnt==0] = np.nan
MHW_mean[MHW_cnt==0] = np.nan
MHW_cum[MHW_cnt==0] = np.nan
MHW_td[MHW_cnt==0] = np.nan
MHW_tc[MHW_cnt==0] = np.nan

# Re-map to run 20E to 380E
#i_20E = np.where(lon_map>20)[0][0]
#lon_map = np.append(lon_map[i_20E:], lon_map[:i_20E]+360)
#SST_mean = np.append(SST_mean[:,i_20E:], SST_mean[:,:i_20E], axis=1)
#MHW_total = np.append(MHW_total[:,i_20E:], MHW_total[:,:i_20E], axis=1)
#MHW_cnt = np.append(MHW_cnt[:,i_20E:], MHW_cnt[:,:i_20E], axis=1)
#MHW_dur = np.append(MHW_dur[:,i_20E:], MHW_dur[:,:i_20E], axis=1)
#MHW_max = np.append(MHW_max[:,i_20E:], MHW_max[:,:i_20E], axis=1)
#MHW_mean = np.append(MHW_mean[:,i_20E:], MHW_mean[:,:i_20E], axis=1)
#MHW_cum = np.append(MHW_cum[:,i_20E:], MHW_cum[:,:i_20E], axis=1)
#MHW_td = np.append(MHW_td[:,i_20E:], MHW_td[:,:i_20E], axis=1)
#MHW_tc = np.append(MHW_tc[:,i_20E:], MHW_tc[:,:i_20E], axis=1)

regions = np.genfromtxt('data/regions', delimiter=',')
regions = np.genfromtxt('data/regions.6oct2017', delimiter=',')
Nreg = regions.shape[0]/5

# Maps

llon, llat = np.meshgrid(lon_map, lat_map)

plt.figure(figsize=(35,4))

plt.clf()
plt.subplot(1,3,3)
proj = bm.Basemap(projection='robin', lon_0=180, resolution='c')
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, MHW_cnt, levels=np.arange(0.5,3.5+0.5,0.5), cmap=plt.cm.YlOrRd)
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
plt.clim(0.75,3.75)
plt.title('Frequency')
h.set_label('annual event count')

#plt.clf()
plt.subplot(1,3,2)
tmp = MHW_dur.copy()
tmp[tmp>40] = 40
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, tmp, levels=[5,10,15,20,30,40], cmap=plt.cm.YlGnBu)
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
h.set_ticklabels(['5', '10', '15', '20', '30', '>40'])
plt.clim(-5,50)
plt.title('Duration')
h.set_label('days')

#plt.clf()
plt.subplot(1,3,1)
proj = bm.Basemap(projection='robin', lon_0=180, resolution='c')
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, MHW_max, levels=[0.,0.5,1,1.5,2,2.5,3,3.5,4,5], cmap=plt.cm.rainbow)
plt.clim(0,4.5)
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
plt.title('Intensity')
h.set_label('$^\circ$C')

# plt.savefig('../../documents/12_Drivers_and_Processes/figures/MHW_int_dur_freq.png', bbox_inches='tight', pad_inches=0.5, dpi=300)

# New map with basic panel

plt.figure(figsize=(15,6))

plt.clf()
plt.subplot(2,2,4)
proj = bm.Basemap(projection='robin', lon_0=180, resolution='c')
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, MHW_cnt, levels=np.arange(0.5,3.5+0.5,0.5), cmap=plt.cm.YlOrRd)
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
plt.clim(0.75,3.75)
plt.title('Frequency')
h.set_label('annual event count')

#plt.clf()
plt.subplot(2,2,3)
tmp = MHW_dur.copy()
tmp[tmp>40] = 40
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, tmp, levels=[5,10,15,20,30,40], cmap=plt.cm.YlGnBu)
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
h.set_ticklabels(['5', '10', '15', '20', '30', '>40'])
plt.clim(-5,50)
plt.title('Duration')
h.set_label('days')

#plt.clf()
plt.subplot(2,2,2)
proj = bm.Basemap(projection='robin', lon_0=180, resolution='c')
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(llon, llat)
plt.contourf(lonproj, latproj, MHW_max, levels=[0.,0.5,1,1.5,2,2.5,3,3.5,4,5], cmap=plt.cm.rainbow)
plt.clim(0,4.5)
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
plt.title('Intensity')
h.set_label('$^\circ$C')

plt.subplot(2,2,1)
proj = bm.Basemap(projection='robin', lon_0=180, resolution='c')
proj.fillcontinents(color='k', lake_color='k')
lonproj, latproj = proj(regions[:,0], regions[:,1])
for i in range(Nreg):
    #plt.plot(lonproj[i*5:((i+1)*5)], latproj[i*5:((i+1)*5)], '-', color='0.4', linewidth=1.5)
    plt.plot(np.append(lonproj[i*5:((i+1)*5)], lonproj[i*5]), np.append(latproj[i*5:((i+1)*5)], latproj[i*5]), '-', color='0.4', linewidth=1.5)
h = plt.colorbar()
plt.title('Regions')
h.set_label('')

# plt.savefig('../../documents/12_Drivers_and_Processes/figures/MHW_int_dur_freq_reg_orig.png', bbox_inches='tight', pad_inches=0.5, dpi=300)


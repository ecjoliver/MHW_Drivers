'''

  Software which uses the MHW definition
  of Hobday et al. (2015) applied to 
  select SST time series around the globe

'''

# Load required modules

import numpy as np
from scipy import io
from scipy import linalg
from scipy import stats
from scipy import interpolate
from scipy import signal
from datetime import date
from netCDF4 import Dataset

from matplotlib import pyplot as plt
import mpl_toolkits.basemap as bm

import marineHeatWaves as mhw
import trendSimAR1

import ecoliver as ecj

#
# observations
#

pathroot = '/mnt/erebor/'
pathroot = '/mnt/EREBOR/'
pathroot = '/mnt/tmp/'
#pathroot = '/bs/projects/geology/Oliver/'
#pathroot = '/home/ecoliver/Desktop/'
header = pathroot+'data/sst/noaa_oi_v2/avhrr/'
file0 = header + '1982/avhrr-only-v2.19820101.nc'
t, dates, T, year, month, day, doy = ecj.timevector([1982, 1, 1], [2017, 12, 31])

#
# lat and lons of obs
#

fileobj = Dataset(file0, mode='r')
lon = fileobj.variables['lon'][:].astype(float)
lat = fileobj.variables['lat'][:].astype(float)
fill_value = fileobj.variables['sst']._FillValue.astype(float)
scale = fileobj.variables['sst'].scale_factor.astype(float)
offset = fileobj.variables['sst'].add_offset.astype(float)
fileobj.close()

#
# Size of mhwBlock variable
#

matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + str(300).zfill(4) + '.mat')
mhws, clim = mhw.detect(t, matobj['sst_ts'][300,:])
mhwBlock = mhw.blockAverage(t, mhws)
years = mhwBlock['years_centre']
NB = len(years)

#
# initialize some variables
#

pctile = 98 # 90 # Percentile for calculation of MHWs
alpha = 0.05
X = len(lon)
Y = len(lat)
#i_which = range(0,X,10)
#j_which = range(0,Y,10)
#i_which = range(0,X,4)
#j_which = range(0,Y,4)
i_which = range(0,X)
j_which = range(0,Y)
DIM = (len(j_which), len(i_which))
SST_mean = np.NaN*np.zeros(DIM)
MHW_total = np.NaN*np.zeros(DIM)
MHW_cnt = np.NaN*np.zeros(DIM)
MHW_dur = np.NaN*np.zeros(DIM)
MHW_max = np.NaN*np.zeros(DIM)
MHW_mean = np.NaN*np.zeros(DIM)
MHW_cum = np.NaN*np.zeros(DIM)
MHW_var = np.NaN*np.zeros(DIM)
MHW_td = np.NaN*np.zeros(DIM)
MHW_tc = np.NaN*np.zeros(DIM)
SST_tr = np.NaN*np.zeros(DIM)
MHW_cnt_tr = np.NaN*np.zeros(DIM)
MHW_dur_tr = np.NaN*np.zeros(DIM)
MHW_max_tr = np.NaN*np.zeros(DIM)
MHW_mean_tr = np.NaN*np.zeros(DIM)
MHW_cum_tr = np.NaN*np.zeros(DIM)
MHW_var_tr = np.NaN*np.zeros(DIM)
MHW_td_tr = np.NaN*np.zeros(DIM)
MHW_tc_tr = np.NaN*np.zeros(DIM)
DIM2 = (len(j_which), len(i_which), 2)
SST_dtr = np.NaN*np.zeros(DIM2)
MHW_cnt_dtr = np.NaN*np.zeros(DIM2)
MHW_dur_dtr = np.NaN*np.zeros(DIM2)
MHW_max_dtr = np.NaN*np.zeros(DIM2)
MHW_mean_dtr = np.NaN*np.zeros(DIM2)
MHW_cum_dtr = np.NaN*np.zeros(DIM2)
MHW_var_dtr = np.NaN*np.zeros(DIM2)
MHW_td_dtr = np.NaN*np.zeros(DIM2)
MHW_tc_dtr = np.NaN*np.zeros(DIM2)
N_ts = np.zeros((len(j_which), len(i_which), NB))
SST_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_cnt_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_dur_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_max_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_mean_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_cum_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_var_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_td_ts = np.zeros((len(j_which), len(i_which), NB))
MHW_tc_ts = np.zeros((len(j_which), len(i_which), NB))
lon_map =  np.NaN*np.zeros(len(i_which))
lat_map =  np.NaN*np.zeros(len(j_which))

# Perform AR1 fit and "excess trends" calculations?
sw_arfit = False
Nens = 500
ar1_tau = np.NaN*np.zeros(DIM)
ar1_sig_eps = np.NaN*np.zeros(DIM)
ar1_putrend_cnt = np.NaN*np.zeros(DIM)
ar1_putrend_mean = np.NaN*np.zeros(DIM)
ar1_putrend_max = np.NaN*np.zeros(DIM)
ar1_putrend_dur = np.NaN*np.zeros(DIM)
ar1_pltrend_cnt = np.NaN*np.zeros(DIM)
ar1_pltrend_mean = np.NaN*np.zeros(DIM)
ar1_pltrend_max = np.NaN*np.zeros(DIM)
ar1_pltrend_dur = np.NaN*np.zeros(DIM)
ar1_mean_cnt = np.NaN*np.zeros(DIM)
ar1_mean_mean = np.NaN*np.zeros(DIM)
ar1_mean_max = np.NaN*np.zeros(DIM)
ar1_mean_dur = np.NaN*np.zeros(DIM)

# Load climate mode indices
sw_noENSO = False
index = {}
# ENSO
dat = np.genfromtxt('../../data/modes/mei.dat', skip_header=10)[:,1:].flatten()
t_mth, tmp, tmp, year_mth, month_mth, day_mth, tmp = ecj.timevector([1950, 1, 1], [2016, 12, 31])
t_mth = t_mth[day_mth==1]
f = interpolate.interp1d(t_mth, dat, fill_value="extrapolate")
index['MEI'] = f(t)
index['MEI'] = index['MEI'] - index['MEI'].mean()
# Regression predictor, Mode index with monthly lead-lags up to a year
# This should be properly called a "Distributed lag [lead-lag] model"
Nl = 12 # Number of lead-lags
Ll = 30 # Length of leads/lags [days]
MPPI = {} # Moore-Penrose pseudo-inverse of X
for mode in index.keys():
    # Old method
    #x = np.array([np.ones(T), index[mode]]).T
    #MPPI[mode] = np.dot(np.linalg.inv(np.dot(x.T, x)), x.T)
    # New method (multiple lead-lags)
    x = np.array([np.ones(T),])
    # Mode leads SST
    for k in range(1,Nl+1):
        tmp = np.append(index[mode][np.arange(T)[k*Ll:]], np.zeros(k*Ll))
        x = np.append(x, np.array([tmp,]), axis=0)
    # 0-lag
    x = np.append(x, np.array([index[mode],]), axis=0)
    # Mode lags SST
    for k in range(1,Nl+1):
        tmp = np.append(np.zeros(k*Ll), index[mode][np.arange(T)[:-k*Ll]])
        x = np.append(x, np.array([tmp,]), axis=0)
    # Hilbert transform method
    #x = np.array([np.ones(T), index[mode], np.imag(signal.hilbert(index[mode]))])

# Theil-Sen trend function
def meanTrend_TS(mhwBlock, alpha=0.05):
    # Initialize mean and trend dictionaries
    mean = {}
    trend = {}
    dtrend = {}
#
    # Construct matrix of predictors, first column is all ones to estimate the mean,
    # second column is the time vector, equal to zero at mid-point.
    t = mhwBlock['years_centre']
    X = t-t.mean()
#
    # Loop over all keys in mhwBlock
    for key in mhwBlock.keys():
        # Skip time-vector keys of mhwBlock
        if (key == 'years_centre') + (key == 'years_end') + (key == 'years_start'):
            continue
#
        # Predictand (MHW property of interest)
        y = mhwBlock[key]
        valid = ~np.isnan(y) # non-NaN indices
#
        # Perform linear regression over valid indices
        if np.sum(~np.isnan(y)) > 0: # If at least one non-NaN value
            slope, y0, beta_lr, beta_up = stats.mstats.theilslopes(y[valid], X[valid], alpha=1-alpha)
            beta = np.array([y0, slope])
        else:
            beta_lr, beta_up = [np.nan, np.nan]
            beta = [np.nan, np.nan]
#
        # Insert regression coefficients into mean and trend dictionaries
        mean[key] = beta[0]
        trend[key] = beta[1]
#
        dtrend[key] = [beta_lr, beta_up]
#
    return mean, trend, dtrend


#
# loop through locations
#

# Tropical Pacific: i = 800, j = 360
# WA: i = 450, j = 260
icnt = 0
for i in i_which:
    print i, 'of', len(lon)-1
#   load SST
    matobj = io.loadmat(header + 'timeseries/avhrr-only-v2.ts.' + str(i+1).zfill(4) + '.mat')
    sst_ts = matobj['sst_ts']
    lon_map[icnt] = lon[i]
#   loop over j
    jcnt = 0
    for j in j_which:
        lat_map[jcnt] = lat[j]
        sst = sst_ts[j,:].copy()
        #if np.logical_not(np.isfinite(sst.sum())): # check for land
        if np.logical_not(np.isfinite(sst.sum())) + ((sst<-1).sum()>0): # check for land, ice
            jcnt += 1
            continue
        # Mode-component of SST
        if sw_noENSO:
            sst_wENSO = sst.copy()
            # Old method
            #beta = np.dot(MPPI['MEI'], np.array([sst,]).T)[1]
            #sst = sst - beta*index['MEI']
            #plt.plot(dates, sst)
            #plt.plot(dates, sst - beta*index['MEI'])
            # New method
            #beta = linalg.lstsq(x.T, sst)[0][1:]
            #plt.plot(dates, sst - np.dot(np.array([beta,]), x[1:])[0,:])
            #sst = sst - np.dot(np.array([beta,]), x[1:])[0,:]
            # Old method on ssta
            #mhws, clim = mhw.detect(t, sst, pctile=pctile)
            #ssta = sst - clim['seas']
            #beta = np.dot(MPPI['MEI'], np.array([ssta,]).T)[1]
            #sst = clim['seas'] + ssta - beta*index['MEI']
            # New method on ssta
            mhws, clim = mhw.detect(t, sst, pctile=pctile)
            ssta = sst - clim['seas']
            beta = linalg.lstsq(x.T, ssta)[0][1:]
            sst = clim['seas'] + ssta - np.dot(np.array([beta,]), x[1:])[0,:]
        # MHW detection
            mhws, clim = mhw.detect(t, sst, pctile=pctile, alternateClimatology=[t, sst_wENSO])
        else:
            mhws, clim = mhw.detect(t, sst, pctile=pctile)
        mhwBlock = mhw.blockAverage(t, mhws, temp=sst)
        # Total count
        MHW_total[jcnt,icnt] = mhwBlock['count'].sum()
        # Mean and trend
        mean, trend, dtrend = mhw.meanTrend(mhwBlock)
        # Mean and trend
        MHW_cnt[jcnt,icnt], MHW_cnt_tr[jcnt,icnt], MHW_cnt_dtr[jcnt,icnt,:] = mean['count'], trend['count'], dtrend['count']
        MHW_dur[jcnt,icnt], MHW_dur_tr[jcnt,icnt], MHW_dur_dtr[jcnt,icnt,:] = mean['duration'], trend['duration'], dtrend['duration']
        MHW_max[jcnt,icnt], MHW_max_tr[jcnt,icnt], MHW_max_dtr[jcnt,icnt,:] = mean['intensity_max_max'], trend['intensity_max_max'], dtrend['intensity_max_max']
        MHW_mean[jcnt,icnt], MHW_mean_tr[jcnt,icnt], MHW_mean_dtr[jcnt,icnt,:] = mean['intensity_mean'], trend['intensity_mean'], dtrend['intensity_mean']
        MHW_cum[jcnt,icnt], MHW_cum_tr[jcnt,icnt], MHW_cum_dtr[jcnt,icnt,:] = mean['intensity_cumulative'], trend['intensity_cumulative'], dtrend['intensity_cumulative']
        MHW_var[jcnt,icnt], MHW_var_tr[jcnt,icnt], MHW_var_dtr[jcnt,icnt,:] = mean['intensity_var'], trend['intensity_var'], dtrend['intensity_var']
        MHW_td[jcnt,icnt], MHW_td_tr[jcnt,icnt], MHW_td_dtr[jcnt,icnt,:] = mean['total_days'], trend['total_days'], dtrend['total_days']
        MHW_tc[jcnt,icnt], MHW_tc_tr[jcnt,icnt], MHW_tc_dtr[jcnt,icnt,:] = mean['total_icum'], trend['total_icum'], dtrend['total_icum']
        SST_mean[jcnt,icnt], SST_tr[jcnt,icnt], SST_dtr[jcnt,icnt,:] = mean['temp_mean'], trend['temp_mean'], dtrend['temp_mean']
        # Time series
        MHW_cnt_ts[jcnt,icnt,:] += mhwBlock['count']
        MHW_dur_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['duration']))[0]] = mhwBlock['duration'][np.where(~np.isnan(mhwBlock['duration']))[0]]
        MHW_max_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_max_max']))[0]] = mhwBlock['intensity_max_max'][np.where(~np.isnan(mhwBlock['intensity_max_max']))[0]]
        MHW_mean_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_mean']))[0]] = mhwBlock['intensity_mean'][np.where(~np.isnan(mhwBlock['intensity_mean']))[0]]
        MHW_cum_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_cumulative']))[0]] = mhwBlock['intensity_cumulative'][np.where(~np.isnan(mhwBlock['intensity_cumulative']))[0]]
        MHW_var_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['intensity_var']))[0]] = mhwBlock['intensity_var'][np.where(~np.isnan(mhwBlock['intensity_var']))[0]]
        MHW_td_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['total_days']))[0]] = mhwBlock['total_days'][np.where(~np.isnan(mhwBlock['total_days']))[0]]
        MHW_tc_ts[jcnt,icnt,np.where(~np.isnan(mhwBlock['total_icum']))[0]] = mhwBlock['total_icum'][np.where(~np.isnan(mhwBlock['total_icum']))[0]]
        N_ts[jcnt,icnt,:] += (~np.isnan(mhwBlock['duration'])).astype(int)
        SST_ts[jcnt,icnt,:] = mhwBlock['temp_mean']
        # Check for excess trends
        if sw_arfit:
            ar1_tau[jcnt,icnt], ar1_sig_eps[jcnt,icnt], trends, means = trendSimAR1.simulate(t, sst, clim['seas'], SST_tr[jcnt,icnt]*10., Nens)
            #ar1_ptrend_cnt[jcnt,icnt] = 1.*(trends['count'] > trend['count']).sum()/len(trends['count'])
            #ar1_ptrend_mean[jcnt,icnt] = 1.*(trends['intensity_mean'] > trend['intensity_mean']).sum()/len(trends['intensity_mean'])
            #ar1_ptrend_dur[jcnt,icnt] = 1.*(trends['duration'] > trend['duration']).sum()/len(trends['duration'])
            ar1_pltrend_cnt[jcnt,icnt] = np.percentile(trends['count'], 2.5)
            ar1_pltrend_mean[jcnt,icnt] = np.percentile(trends['intensity_mean'], 2.5)
            ar1_pltrend_max[jcnt,icnt] = np.percentile(trends['intensity_max_max'], 2.5)
            ar1_pltrend_dur[jcnt,icnt] = np.percentile(trends['duration'], 2.5)
            ar1_putrend_cnt[jcnt,icnt] = np.percentile(trends['count'], 97.5)
            ar1_putrend_mean[jcnt,icnt] = np.percentile(trends['intensity_mean'], 97.5)
            ar1_putrend_max[jcnt,icnt] = np.percentile(trends['intensity_max_max'], 97.5)
            ar1_putrend_dur[jcnt,icnt] = np.percentile(trends['duration'], 97.5)
            ar1_mean_cnt[jcnt,icnt] = means['count'].mean()
            ar1_mean_mean[jcnt,icnt] = means['intensity_mean'].mean()
            ar1_mean_max[jcnt,icnt] = means['intensity_max_max'].mean()
            ar1_mean_dur[jcnt,icnt] = means['duration'].mean()
        # Up counts
        jcnt += 1
    icnt += 1
    # Save data so far
    #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.lores.N1000'
    #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2015'
    #outfile = '/bs/projects/geology/Oliver/data/MHWs/Trends/mhw_census_TS.2015'
    if sw_noENSO:
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores.noENSO.leadLag.1yr'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores.noENSO.lag0'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores.noENSO.Hilbert'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores.noENSO.Hilbert.ssta'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores.noENSO.lag0.ssta'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.noENSO.leadLag.1yr.monthly.ssta'
        outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores.noENSO.leadLag.1yr.monthly.ssta'
    elif sw_arfit:
        outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.excessTrends.lores'
    else:
        outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2017'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.p98'
        #outfile = '/home/ecoliver/Desktop/data/MHWs/Trends/mhw_census.2016.lores'
    if (i % 100) + (i == i_which[-1]):
        np.savez(outfile, lon_map=lon_map, lat_map=lat_map, SST_mean=SST_mean, MHW_total=MHW_total, MHW_cnt=MHW_cnt, MHW_dur=MHW_dur, MHW_max=MHW_max, MHW_mean=MHW_mean, MHW_cum=MHW_cum, MHW_var=MHW_var, MHW_td=MHW_td, MHW_tc=MHW_tc, SST_tr=SST_tr, MHW_cnt_tr=MHW_cnt_tr, MHW_dur_tr=MHW_dur_tr, MHW_max_tr=MHW_max_tr, MHW_mean_tr=MHW_mean_tr, MHW_cum_tr=MHW_cum_tr, MHW_var_tr=MHW_var_tr, MHW_td_tr=MHW_td_tr, MHW_tc_tr=MHW_tc_tr, SST_dtr=SST_dtr, MHW_cnt_dtr=MHW_cnt_dtr, MHW_dur_dtr=MHW_dur_dtr, MHW_max_dtr=MHW_max_dtr, MHW_mean_dtr=MHW_mean_dtr, MHW_cum_dtr=MHW_cum_dtr, MHW_var_dtr=MHW_var_dtr, MHW_td_dtr=MHW_td_dtr, MHW_tc_dtr=MHW_tc_dtr, SST_ts=SST_ts, MHW_cnt_ts=MHW_cnt_ts, MHW_dur_ts=MHW_dur_ts, MHW_max_ts=MHW_max_ts, MHW_mean_ts=MHW_mean_ts, MHW_cum_ts=MHW_cum_ts, MHW_var_ts=MHW_var_ts, MHW_td_ts=MHW_td_ts, MHW_tc_ts=MHW_tc_ts, N_ts=N_ts, years=years, alpha=alpha, ar1_tau=ar1_tau, ar1_sig_eps=ar1_sig_eps, ar1_putrend_cnt=ar1_putrend_cnt, ar1_putrend_mean=ar1_putrend_mean, ar1_putrend_max=ar1_putrend_max, ar1_putrend_dur=ar1_putrend_dur, ar1_pltrend_cnt=ar1_pltrend_cnt, ar1_pltrend_mean=ar1_pltrend_mean, ar1_pltrend_max=ar1_pltrend_max, ar1_pltrend_dur=ar1_pltrend_dur, ar1_mean_cnt=ar1_mean_cnt, ar1_mean_mean=ar1_mean_mean, ar1_mean_max=ar1_mean_max, ar1_mean_dur=ar1_mean_dur)



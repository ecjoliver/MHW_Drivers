from netCDF4 import Dataset
import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.basemap import Basemap
from datetime import date
import os
import shutil
import datetime
import sys

# NB only works from the command line e.g. python regions_MHW_pc90.py regional_data/avhrr-only-v2.0to30.-10to10.1981.2016.nc
# as need to pass file as argument
# Files will be run on Raijin at /g/data1a/e14/asg561/MHW

# Load marineHeatWaves definition module 
import marineHeatWaves90 as mmhw

# load input file name from parallel_process_MHW_regional_blocks.sh (batch calling routine) 
for arg in sys.argv:
 file=arg

os.chdir('/g/data1a/e14/asg561/MHW')
# load file data

fh = Dataset(file, mode='r')
lon = fh.variables['lon'][:]
lat = fh.variables['lat'][:]
time = fh.variables['time'][:]
t=date(1978,1,1).toordinal()+time
t=t.astype(int)


sst_units = fh.variables['sst'].units
SEVERITY=np.zeros((len(t),len(lat)/8,len(lon)/8)) # severity 
ANOMALY=np.zeros((len(t),len(lat)/8,len(lon)/8))  # SSTA 
CLIMATOLOGY=np.zeros((366,len(lat)/8,len(lon)/8))  # long term mean daily climatology (smoothed)
CLIMATOLOGYpc90=np.zeros((366,len(lat)/8,len(lon)/8))  # 90th percentile climatology (smoothed)
LON=np.zeros((len(lon)/8))
LAT=np.zeros((len(lat)/8))

ANOMALY[:] = np.NAN
SEVERITY[:] = np.NAN
CLIMATOLOGY[:] = np.NAN
CLIMATOLOGYpc90[:] = np.NAN

# destination filename
fstart, fend=file.split(".",1)


dst='/g/data1a/e14/asg561/MHW/mhw_data_90pc_2degree/mhw_severity.pc90.'+fend
ti = np.where((t >= date(1984,1,1).toordinal()) & (t <= date(1984,12,31).toordinal()) )[0]
if not os.path.isfile(dst):
    print(dst) 
    #for mylat in range(0,5):
    for mylat in range(0,len(lat)/8):    
        #sst_lat = fh.variables['sst'][ : , 0, mylat,: ]
        # convert masked array to normal array with filled values
        #if type(sst_lat).__name__ == 'MaskedArray':
        #    sst_lat=sst_lat.filled()

        # remove data prior to 1982
        #sst_lat=sst_lat[ti,:]
        tmp=lat[mylat*8:mylat*8+8]
        LAT[mylat]=tmp.mean()
        for mylon in range(0,len(lon)/8):
            tmp=lon[mylon*8:mylon*8+8]
            LON[mylon]=tmp.mean()
            sst_8 = fh.variables['sst'][ : , 0, mylat*8:mylat*8+8,mylon*8:mylon*8+8 ]
            #if type(sst_8).__name__ == 'MaskedArray':
            #	sst_8=sst_8.filled()
            sst=sst_8.mean(1)
            sst=sst.mean(1)
            # check if sst is empty
            #if not all(x == sst[0] for x in sst):
            if  not isinstance(sst,np.ma.MaskedArray):
                # detrend sst but retain mean
                tiY = np.where((t >= date(1982,1,1).toordinal()) & (t < date(2016,1,1).toordinal()) )[0] # only inlude complete years to calculate trend
                model = np.polyfit(t[tiY], sst[tiY], 1)
                predicted = np.polyval(model, t)
                sst=sst-predicted+sst.mean(0)
                mhws, clim = mmhw.detect(t, sst)
                anomaly=sst-clim['seas'] # ssta 

                # creat time series MHWbool that is 1 for MHW, NaN otherwise
                N = mhws['n_events']
                MHWbool = np.zeros((len(t)))
                MHWbool[:] = np.nan
                for ev0 in np.arange(0, N, 1):
                        t1 = mhws['index_start'][ev0]
                        t2 = mhws['index_end'][ev0]
                        MHWbool[t1:t2+1] = 1

                sev=clim['thresh']-clim['seas']
                sev=np.divide(anomaly,sev)
                # only retain MHW events
                sev = sev * MHWbool


                #sev=clim['thresh']-clim['seas']
                #sev=np.divide(anomaly,sev)
                #sev[sev<1] = np.nan #remove non MHW regions 
                SEVERITY[:,mylat,mylon]= sev
                ANOMALY[:,mylat,mylon]= anomaly
                # isolate leap year for saving
                #ti = np.where((t >= date(1984,1,1).toordinal()) & (t <= date(1984,12,31).toordinal()) )[0]
                tmp=clim['seas']
                CLIMATOLOGY[:,mylat,mylon]= tmp[ti] # extract 1 (leap) year
                tmp=clim['thresh']
                CLIMATOLOGYpc90[:,mylat,mylon]= tmp[ti] # extract 1 (leap) year
                del tmp, clim, mhws, anomaly
            ## END mylon loop
    ## END mylat loop 
    #del tmp, clim, mhws, anomaly 
    

    # open a new netCDF file for writing.
    ncfile = Dataset(dst,'w')
    # create the x and y dimensions.
    ncfile.history = 'Created in regional_MHW_pc90_reducedFileSize_2degree.py pctile=90'+str(datetime.datetime.now())
    ncfile.createDimension('lon',len(LON))
    ncfile.createDimension('lat',len(LAT))
    ncfile.createDimension('time',None)
    ncfile.createDimension('time366',None)

    # create the variable .
    severity = ncfile.createVariable('severity','f4',('time','lat','lon'),zlib=True, least_significant_digit=3)
    ssta = ncfile.createVariable('ssta','f4',('time','lat','lon'),zlib=True, least_significant_digit=4)
#    severity = ncfile.createVariable('severity','f4',('time','lat','lon'))
#    ssta = ncfile.createVariable('ssta','f4',('time','lat','lon'))
    climatology90 = ncfile.createVariable('climatology90','f4',('time366','lat','lon'))
    climatology = ncfile.createVariable('climatology','f4',('time366','lat','lon'))
    lats = ncfile.createVariable('lat','f4',('lat',))
    lons = ncfile.createVariable('lon','f4',('lon',))
    time = ncfile.createVariable('time','f4',('time',))
    time366 = ncfile.createVariable('time366','f4',('time366',))
    lats.units = 'degrees_north'
    lons.units = 'degrees_east'
    time.units = 'days since 01-01-01 00:00:00'
    time366.units = 'days since 01-01-01 00:00:00' # need to correct this
    ssta.units='degrees C'
    ssta.long_name='Daily sea surface temperature anomaly'
    climatology90.units='degrees C'
    climatology90.long_name='90pc Climatological sea surface temperature'
    climatology.units='degrees C'
    climatology.long_name='Climatological sea surface temperature'
    severity.units='ratio SSTA/(CLIM90-CLIMseas)'
    severity.long_name='MHW Severity'

    # write data to variable.
    lats[:] = LAT
    lons[:] = LON
    time[:] = t-1 # for some reason dates are out by 1 day; minus 1 to correct
    time366[:] = t[ti]-1  # for some reason dates are out by 1 day; minus 1 to correct
    severity[:] = SEVERITY
    ssta[:] = ANOMALY
    climatology90[:] = CLIMATOLOGYpc90
    climatology[:] = CLIMATOLOGY
    # close the file.
    ncfile.close()
    print('*** SUCCESS writing netcdf file ***')
    ## END if dst file exists
    del ssta,  ANOMALY, ncfile, fh


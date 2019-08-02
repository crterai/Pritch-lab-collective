import numpy as np
import scipy.io as sio
import random as rand
import multiprocessing
import time
import pickle

from multiprocessing import Pool, Process
from joblib import Parallel, delayed
from scipy.stats import genextreme as gev
from netCDF4 import Dataset

tic = time.time()
# Import CaMa grid data 
CaMaFile = Dataset('/work/04268/tg835671/stampede2/PythonBootstrap/fldare1161.nc')
lat =  CaMaFile.variables['lat'][:] 
lon = CaMaFile.variables['lon'][:] 

#Import matlab arrays of maximum annual outflow 
#matData   = sio.loadmat('/work/04268/tg835671/stampede2/PythonBootstrap/control_outflw.mat')
#maxYearly = matData['maxYearly_full']

#fullData   = sio.loadmat('/work/04268/tg835671/stampede2/PythonBootstrap/full_outflw.mat')
#fullYearly = fullData['maxYearly_full']

#physData = sio.loadmat('/work/04268/tg835671/stampede2/PythonBootstrap/physiology_outflw.mat')
#physYearly = physData['maxYearly_full']

#radData = sio.loadmat('/work/04268/tg835671/stampede2/PythonBootstrap/radiation_outflw.mat')
#radYearly = radData['maxYearly_full']

toc = time.time() - tic
print("Time elapsed reading in data: ", toc)


# Define return period and number of bootstrap samples 
retPeriod = 100.0
nsamples = 1000

# Define empty arrays for storing data in
Tfull_boot = np.full([nsamples],np.nan)
Tphys_boot = np.full([nsamples],np.nan)
Trad_boot  = np.full([nsamples],np.nan)

#Tfull_Mean = np.full([len(lon),len(lat)],np.nan)
#Tfull_2_5pct = np.full([len(lon),len(lat)],np.nan)
#Tfull_5pct   = np.full([len(lon),len(lat)],np.nan)
#Tfull_15pct  = np.full([len(lon),len(lat)],np.nan)
#Tfull_85pct  = np.full([len(lon),len(lat)],np.nan)
#Tfull_95pct  = np.full([len(lon),len(lat)],np.nan)
#Tfull_97_5pct = np.full([len(lon),len(lat)],np.nan)

#Tphys_Mean    = np.full([len(lon),len(lat)],np.nan)
#Tphys_2_5pct  = np.full([len(lon),len(lat)],np.nan)
#Tphys_5pct    = np.full([len(lon),len(lat)],np.nan)
#Tphys_15pct   = np.full([len(lon),len(lat)],np.nan)
#Tphys_85pct   = np.full([len(lon),len(lat)],np.nan)
#Tphys_95pct   = np.full([len(lon),len(lat)],np.nan)
#Tphys_97_5pct = np.full([len(lon),len(lat)],np.nan)

#Trad_Mean    = np.full([len(lon),len(lat)],np.nan)
#Trad_2_5pct  = np.full([len(lon),len(lat)],np.nan)
#Trad_5pct    = np.full([len(lon),len(lat)],np.nan)
#Trad_15pct   = np.full([len(lon),len(lat)],np.nan)
#Trad_85pct   = np.full([len(lon),len(lat)],np.nan)
#Trad_95pct   = np.full([len(lon),len(lat)],np.nan)
#Trad_97_5pct = np.full([len(lon),len(lat)],np.nan)

def computeGEV(recArray,recFullArr,recPhysArr,recRadArr):
    iRand = np.random.randint(0,30,size=30)
            
    rec     = recArray[iRand]     # Sub-sample by random integers
    recFull = recFullArr[iRand]
    recPhys = recPhysArr[iRand]
    recRad  = recRadArr[iRand]
       
    # Check for NaN 
    validCount = np.count_nonzero(~np.isnan(rec))
          
    if validCount>5: 
        #print("Starting to compute...")
            
        gev_fit = gev.fit(rec)
        K = gev_fit[0]
        M = gev_fit[1]
        S = gev_fit[2]
    
        #Compute discharge of 100-year flood 
        R100 = gev.ppf((1-(1/retPeriod)),K,loc=M,scale=S)
    
        #Compute fit for CO2 experiments
        gevFull = gev.fit(recFull)
        gevPhys = gev.fit(recPhys)
        gevRad  = gev.fit(recRad)
    
        cdfFull = gev.cdf(R100,*gevFull)    #Compute CDF at the point of R100 
        cdfPhys = gev.cdf(R100,*gevPhys)
        cdfRad  = gev.cdf(R100,*gevRad)
    
        TfullInt = 1/(1-cdfFull)    #Return Period 
        TphysInt = 1/(1-cdfPhys)
        TradInt  = 1/(1-cdfRad)

        #print("Finished with a bootstrap...")

    else: 
        TfullInt = np.nan
        TphysInt = np.nan
        TradInt = np.nan 
        #print("Finished computing GEV...")    

    return TfullInt, TphysInt, TradInt

import os
lonInd = int(os.getenv('XINDEX_CURRENT'))
latInd = int(os.getenv('YINDEX_CURRENT'))

print("lonInd: ", lonInd)
print("latInd: ", latInd)

for i in [lonInd]:
    for j in [latInd]: 
#for i in range(len(lon)/7):
#    for j in range(len(lat)/7):

        # READ in location-specific data
        maxYearlyFile = '/work/04268/tg835671/stampede2/PythonBootstrap/maxYearlyFiles/maxYearly_lonIndex'+str(i)+'_latIndex'+str(j)+'.pkl' 
        [recO,recFullO,recPhysO,recRadO]= pickle.load( open( maxYearlyFile, "rb" ) )
        
        #recO     = maxYearly[i,j,:]   # Read in original data arrays at gridcell 
        #recFullO = fullYearly[i,j,:]
        #recPhysO = physYearly[i,j,:]
        #recRadO  = radYearly[i,j,:]
 
        print("Longitude index being used: ", i)

        # Check for NaN 
        validCount = np.count_nonzero(~np.isnan(recO))

        if validCount<5:
            Tfull_boot[:] = np.nan
            Tphys_boot[:] = np.nan
            Trad_boot[:]  = np.nan
            
            Tfull_Mean   = np.nan
            Tfull_2_5pct = np.nan
            Tfull_5pct   = np.nan
            Tfull_15pct  = np.nan
            Tfull_85pct  = np.nan
            Tfull_95pct  = np.nan
            Tfull_97_5pct = np.nan

            Tphys_Mean    = np.nan
            Tphys_2_5pct  = np.nan
            Tphys_5pct    = np.nan
            Tphys_15pct   = np.nan
            Tphys_85pct   = np.nan
            Tphys_95pct   = np.nan
            Tphys_97_5pct = np.nan

            Trad_Mean    = np.nan
            Trad_2_5pct  = np.nan
            Trad_5pct    = np.nan
            Trad_15pct   = np.nan
            Trad_85pct   = np.nan
            Trad_95pct  = np.nan
            Trad_97_5pct = np.nan

            print("No valid points for lon index ",i, "and lat index ", j)
        else: 
            # Compute nsamples-member boostrap array of GEV return periods for single gridcell using 10 workers via Parallel
            tic = time.time()
            print("Beginning to compute bootstraps... ")        
            for b in range(nsamples):
                [full,phys,rad] = computeGEV(recO,recFullO,recPhysO,recRadO) 
                Tfull_boot[b] = full
                Tphys_boot[b] = phys
                Trad_boot[b]  = rad
            #Tbig = Parallel(n_jobs=10)(delayed(computeGEV)(recO,recFullO,recPhysO,recRadO) for b in range(nsamples))       
            toc = time.time() - tic
        
            if Tfull_boot[0]>=0:    #Print timing stats to screen 
                print("Time elapsed in bootstrap GEV: ", toc)

            #Place returned tuples in appropriate arrays 
            #for b in range(len(Tbig)):
            #   Tfull_boot[b] = Tbig[b][0]
            #    Tphys_boot[b] = Tbig[b][1]
            #    Trad_boot[b]  = Tbig[b][2]
        
            # Compute mean and perentiles of bootstraps
            Tfull_Mean   = np.nanmean(Tfull_boot,axis=0)
            Tfull_2_5pct = np.nanpercentile(Tfull_boot,2.5,axis=0)
            Tfull_5pct   = np.nanpercentile(Tfull_boot,5,axis=0)
            Tfull_15pct  = np.nanpercentile(Tfull_boot,15,axis=0)
            Tfull_85pct  = np.nanpercentile(Tfull_boot,85,axis=0)
            Tfull_95pct  = np.nanpercentile(Tfull_boot,95,axis=0)
            Tfull_97_5pct = np.nanpercentile(Tfull_boot,97.5,axis=0)

            Tphys_Mean    = np.nanmean(Tphys_boot,axis=0)
            Tphys_2_5pct  = np.nanpercentile(Tphys_boot,2.5,axis=0)
            Tphys_5pct    = np.nanpercentile(Tphys_boot,5,axis=0)
            Tphys_15pct   = np.nanpercentile(Tphys_boot,15,axis=0)
            Tphys_85pct   = np.nanpercentile(Tphys_boot,85,axis=0)
            Tphys_95pct   = np.nanpercentile(Tphys_boot,95,axis=0)
            Tphys_97_5pct = np.nanpercentile(Tphys_boot,97.5,axis=0)

            Trad_Mean    = np.nanmean(Trad_boot,axis=0)
            Trad_2_5pct  = np.nanpercentile(Trad_boot,2.5,axis=0)
            Trad_5pct    = np.nanpercentile(Trad_boot,5,axis=0)
            Trad_15pct   = np.nanpercentile(Trad_boot,15,axis=0)
            Trad_85pct   = np.nanpercentile(Trad_boot,85,axis=0)
            Trad_95pct   = np.nanpercentile(Trad_boot,95,axis=0)
            Trad_97_5pct = np.nanpercentile(Trad_boot,97.5,axis=0)

        # WRITE DATA TO FILE 
        #tic1 =time.time()
        #Save as pickle 
        fileName = '/work/04268/tg835671/stampede2/PythonBootstrap/bootstrapResults-NoParallel/testBoot1000_lonIndex'+str(i)+'_latIndex'+str(j)+'.pkl' 
        f_myFile = open(fileName,'wb')
        pickle.dump([Tfull_Mean, Tfull_2_5pct, Tfull_5pct, Tfull_15pct, Tfull_85pct, Tfull_95pct, Tfull_97_5pct, Tphys_Mean, Tphys_2_5pct, Tphys_5pct, Tphys_15pct, Tphys_85pct, Tphys_95pct, Tphys_97_5pct, Trad_Mean, Trad_2_5pct, Trad_5pct, Trad_15pct, Trad_85pct, Trad_95pct, Trad_97_5pct], f_myFile)
        f_myFile.close()  
        #toc1=time.time() - tic1
        print("File saved for lon ",i, " and lat ",j)


#import os
#nboot = os.getenv('iboot')
#fileName = '/gdata/pritchard2/mdfowler/Flooding-physiology/Python/results_testBootstrap/pyTboot_' + str(nboot) + '.mat'


 


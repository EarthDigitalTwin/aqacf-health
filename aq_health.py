#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 14:23:54 2022

@author: schakrab
"""

import numpy as np
import glob
import datetime
from datetime import date
from datetime import timedelta
from dateutil.relativedelta import relativedelta
from datetime import datetime
import pandas as pd
import os
from netCDF4 import Dataset
from osgeo import gdal
import geopandas as gpd
from shapely.geometry import MultiPolygon, Polygon, Point
from scipy import io
# from pm_health_calc import pm_health_calc
# from no_health_calc import no_health_calc
# from o3_health_calc import o3_health_calc

def inpolygon(cc,dd,xx,yy,multiplgn):

   temp1=np.zeros((cc,dd))
   for kk in range(0,1):
       for ii in range(cc):
            for jj in range(dd):
                if Point((xx[ii,jj],yy[ii,jj])).within(multiplgn) and temp1[ii,jj] == 0:
                 # print('inside',(xx[ii,jj],yy[ii,jj])) 
                  temp1[ii,jj]=1 
                  #print(np.sum(temp1))
   return(temp1)


# Path to script
base_dir = '/Users/schakrab/Desktop/Paper/input_folder/'
#City Initials
INIT     = 'HOU'
#City Name
NAME     = 'Houston'
#State city is in
STATE    = 'Texas'
SEX = 'both'
#Set-Up
########################################################################################################################################################
#command lines### list of names of the cities and state ### Hardparse
########################################################################################################################################################
#Directories

out_dir  = base_dir+'outputs/'
in_dir   = base_dir+'static/'
sens_dir = base_dir+'inputs/'+INIT+'/'
emis_dir = in_dir+'emissions/'
othe_dir = in_dir+'other/'
ozon_dir = in_dir+'ozone/'
nei_dir  = in_dir+'nei/'
hlth_dir = in_dir+'health/'



#Files( state if teh file is missing or corrupted)
# make a structure of teh imput files.
pmfile     = othe_dir+'V5GL02.HybridPM25.Global.201101-201112.nc' #monthly PM data
nofile     = othe_dir+'TROP_2011_DC_US_v2.nc' # NO data
popfile    = othe_dir+'gpw_v4_population_count_rev11_2010_30_sec.tif'# population data
maskfile   = othe_dir+'cities_13k_largerext.tif' # mask data
lookupfile = othe_dir+'Revised_GHS-SMOD HDC_Lookup_DAM_112421.csv' # infor about the cities
areafile   = othe_dir+'GC_M2_US_05x0667.nc' # area of each grid
#fh = Dataset(areafile, mode='r')
gbdpopfile = hlth_dir+'IHME_GBD_2019_POP_SYA_2011_Y2021M01D28.CSV' # global disease burden pop for states
ctypopfile = hlth_dir+'USCITIES_07212022_Data.csv' # read cities pop and kids pop
mortfile   = hlth_dir+'IHME-GBD_2019_DATA.csv' # death data age group wise
asthfile   = hlth_dir+'IHME-GBD_2019_DATA_ASTHMA.csv'# asthma data age group wise
# file containing codes for age+ location+ and cause
acodefile = hlth_dir+'IHME_GBD_2019_CONTEXTS_BY_AGE_Y2020M12D18.csv'
lcodefile = hlth_dir+'IHME_GBD_2019_ALL_LOCATIONS_HIERARCHIES_Y2022M01D20.csv'
ccodefile = hlth_dir+'IHME_GBD_2019_CAUSE_HIERARCHY_Y2020M11D25.csv'
files_o3 = glob.glob(ozon_dir+'/*.txt') #
flen_o3 = len(files_o3)

# Assignments
IIPAR     = 91
JJPAR     = 89
dmax      = 396
k_pm      = 34
k_o3      = 62
k_no      = 34
ks        = [k_pm,k_o3,k_no]

di        = (date.toordinal(date(2010, 11, 30))) #matlab starts from 0000 python from 0001
pollutant = ['PM','O3','NO']
pmax      = len(pollutant)
# emission factors
SPD        = 60 * 60 * 24
CM2M       = 10000
MW         = [46.005, 12.0107+15.9994,12.0107]
KGCONV     = np.array(MW) * 0.001 / (6.0221415E23)
efactor    = SPD * CM2M * KGCONV
# health information
ages       = [25,30,35,40,45,50,55,60,65,70,75,80,85,90,95]
ageids     = [10,11,12,13,14,15,16,17,18,19,20,30,31,32,235]
amax       = len(ages)
# mortality rate outcome names
mortid      = [976,509,426,322,494,493]
mmax        = len(mortid)     


# Initializations
# Sensitivities
PMGDT = np.zeros((IIPAR,JJPAR,dmax,k_pm)) 
PMEMS = PMGDT
O3GDT = np.zeros((IIPAR,JJPAR,dmax,k_o3)) 
O3EMS = O3GDT
NOGDT = np.zeros((IIPAR,JJPAR,dmax,k_no)) 
NOEMS = NOGDT

# Other
inames    = np.zeros((62,3))
O3        = np.zeros((91,89)) 
CT = np.zeros((91,89))

inames0= []  
inames1= [] 
inames2= []

# (1) Read in adjoint sensitivities and aggregate to annual time-scale
print('Begin adjoint sensitivity read and aggregation')

# loop through three pollutants

for p  in range(pmax):
    A   = np.zeros((IIPAR,JJPAR,dmax,ks[p])) #######################################
    G = np.zeros((IIPAR,JJPAR,dmax,ks[p]))#######################################
    TOT = np.zeros((ks[p])) #######################################
    C = 0

    #loop through months    
    for m in range(1,13):#1,13
      # get date quantities
        if m <10:# mon  = pad(str(m),2,'left','0') 
         mon = str(m).rjust(1 + len(str(m)), '0')
        if m>9: 
         mon = str(m)#.rjust(1 + len(str(m)), '0')
       
        dtObj = datetime.strptime(str(date(2011,m,1)), '%Y-%m-%d')
        past_date = dtObj - pd.DateOffset(months=1)
        d0= (date.toordinal(past_date.date()))-di-1
        post_date = dtObj + pd.DateOffset(months=1)- pd.DateOffset(days=1)
        df= (date.toordinal(post_date.date()))-di
        dlen = (df-d0)
        #print(p,m,dtObj,past_date,post_date,d0,df,dlen)
       #get paths to sensitivity files and skip empty months
        adj =base_dir+'inputs/'+INIT+'/'+'ems.adj.'+mon+'.'+INIT+'.'+pollutant[p]+'.nc'
        gdt =base_dir+'inputs/'+INIT+'/'+'gctm.gdt.'+mon+'.'+INIT+'.'+pollutant[p]+'.nc'
        #add file missing 
        if os.path.isfile(adj):
            fhadj = Dataset(adj, mode='r')
            ivar=fhadj.variables.keys() 
            fhgdt = Dataset(gdt, mode='r')
            ivargdt=fhgdt.variables.keys() 
            
            ilen=len(ivar)
            ivar=list(ivar)
           # print(adj,np.shape(ivar))
            if np.any(TOT !=0):
                C=1
            for  i in range(6,ilen) : 
                    if C==1:
                     if TOT[i-6] ==0:
                         continue
                    name=ivar[i]
                    #print(name)  
                    if p==0 :
                     inames0.append(['IJ-GDE-S__'+name[9:len(name)]])
                     tmpg = np.transpose(fhgdt.variables['IJ-GDE-S__'+name[9:len(name)]][:])
                     if np.sum(tmpg) == 0:
                         continue
                     
                     tmpa =  np.transpose(fhadj.variables[name][:])
                     A[:,:,d0:df,i-6]=A[:,:,d0:df,i-6]+tmpa[:,:,0:dlen]*dlen
                     G[:,:,d0:df,i-6] = G[:,:,d0:df,i-6] + tmpg[:,:,0:dlen] 
                     TOT[i-6]=np.sum(tmpg) 
                    if p==1 :
                     inames1.append(['IJ-GDE-S__'+name[9:len(name)]])
                     tmpg = np.transpose(fhgdt.variables['IJ-GDE-S__'+name[9:len(name)]][:])
                     if np.sum(tmpg) == 0:
                         continue
                     tmpa =  np.transpose(fhadj.variables[name][:])
                     A[:,:,d0:df,i-6]=A[:,:,d0:df,i-6]+tmpa[:,:,0:dlen]*dlen
                     G[:,:,d0:df,i-6] = G[:,:,d0:df,i-6] + tmpg[:,:,0:dlen] 
                     TOT[i-6]=np.sum(tmpg)
                    if p==2 :
                     inames2.append(['IJ-GDE-S__'+name[9:len(name)]])
                     tmpg = np.transpose(fhgdt.variables['IJ-GDE-S__'+name[9:len(name)]][:])
                     if np.sum(tmpg) == 0:
                         continue
                     tmpa =  np.transpose(fhadj.variables[name][:])
                     A[:,:,d0:df,i-6]=A[:,:,d0:df,i-6]+tmpa[:,:,0:dlen]*dlen
                     G[:,:,d0:df,i-6] = G[:,:,d0:df,i-6] + tmpg[:,:,0:dlen] 
                     TOT[i-6]=np.sum(tmpg)
        #print(m,p,np.max(A))    
            fhadj.close()
            fhgdt.close()  
            
            # close i loop
    # close month looop     
    #print(np.max(A))
    if p == 0:
     PMEMS = A
     PMGDT = G
    if p == 1:
     O3EMS = A
     O3GDT = G
    if p == 2:
     NOEMS = A
     NOGDT = G
    
#clopse p loop            
inames0=list(inames0[0:34])
inames1=list(inames1[0:62])
inames2=list(inames2[0:34])            

# coordinate and area data from GEOS-Chem
fhadj = Dataset(adj, mode='r')
latgc = fhadj.variables['LAT'][:]
longc = fhadj.variables['LON'][:]
fhadj.close()
fharea = Dataset(areafile, mode='r')
GCM2 = fharea.variables['DXYP__DXYP'][:]
fharea.close()

#%% (2) Get cost-function data
# % Read in mask data
#-------------------------------------------------------------------------%
dataset = gdal.Open(maskfile, gdal.GA_ReadOnly)#gdal.GA_ReadOnly
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    A = band.ReadAsArray()
       
df = pd.read_csv(lookupfile,skiprows=1, names=['ID_HDC_G0','City_Corrected','Cluster_Corrected','Country_Corrected','GBD_Region_Corrected'])
T=np.squeeze([[df.ID_HDC_G0],[df.City_Corrected],[df.Cluster_Corrected],[df.Country_Corrected],[df.GBD_Region_Corrected]])
for ii in range(len(T[2,:])):
    if T[1,ii] == str(NAME):
        ind=ii+1
A[A!=ind] = 0   
latf = np.linspace(75,-56,15720)
lonf = np.linspace(-179.5000,179.0000,43020)
ymin=np.min(np.where(np.sum(A,axis=1)>0))
ymax=np.max(np.where(np.sum(A,axis=1)>0))
xmin=np.min(np.where(np.sum(A,axis=0)>0))
xmax=np.max(np.where(np.sum(A,axis=0)>0))
# limit to just the city boundaries
Z  = A[ymin:ymax+1,xmin:xmax+1]
la = latf[ymin:ymax+1]
lo = lonf[xmin:xmax+1]

#%-------------------------------------------------------------------------%



# Read in O3 data
#-------------------------------------------------------------------------%

for i in range(len(files_o3)):
    with open(files_o3[i]) as f:
     f=[x.strip() for x in f if x.strip()]
     temp=[tuple(map(float,x.split())) for x in f[0:]]
     temp=np.array(temp)
     temp=np.transpose(temp)  
     O3[np.where(temp>0)]=O3[np.where(temp>0)]+temp[np.where(temp>0)] 
     CT[np.where(temp>0)]  =CT[np.where(temp>0)]+1  
O3 = O3 / CT 
# #O3=np.transpose(O3)
#-------------------------------------------------------------------------%


# Read in PM Data
#-------------------------------------------------------------------------%
PM = Dataset(pmfile, mode='r').variables['GWRPM25'][:]
latpm    = Dataset(pmfile, mode='r').variables['lat'][:]
lonpm = Dataset(pmfile, mode='r').variables['lon'][:]
PM=np.transpose(PM)
#-------------------------------------------------------------------------%

# % Read in Pop data
# %-------------------------------------------------------------------------%
dataset = gdal.Open(popfile, gdal.GA_ReadOnly)#gdal.GA_ReadOnly
for x in range(1, dataset.RasterCount + 1):
    band = dataset.GetRasterBand(x)
    POP = band.ReadAsArray()
latpo  = np.linspace(-90,89.999999,21600)
lonpo = np.linspace(-180,180,43200)
POP=np.flipud(POP)
# %-------------------------------------------------------------------------%


# Read in NO2 data
#-------------------------------------------------------------------------%
NO = Dataset(nofile, mode='r').variables['DS'][:]
latno    = Dataset(nofile, mode='r').variables['lat'][:]
lonno= Dataset(nofile, mode='r').variables['lon'][:]
NO=np.transpose(NO)
#-------------------------------------------------------------------------%

#Calculate city mask at 0.5° x 0.667° resolution
#-------------------------------------------------------------------------%
CITY_MASK = np.zeros((91,89))
for i in range(0,len(lo)):
 for j in range(0,len(la)):
    x=((abs(longc-lo[i])))
    x=np.where(x==np.min(x))
    y=((abs(latgc-la[j])))
    y=np.where(y==np.min(y))
    if Z[j,i] != 0: 
          CITY_MASK[x,y] = CITY_MASK[x,y] + 1/4840
CITY_MASK = np.transpose(CITY_MASK)      ########################    
#check value with matlabprint(np.sum(Z,axis=0) )       
#-------------------------------------------------------------------------% 


# Calculate state mask at at 0.5° x 0.667° resolution
#-------------------------------------------------------------------------%
# read in state shape file

S = gpd.read_file(othe_dir+'s_22mr22.shp')
Sgeo=S.geometry; XS=[]; YS=[]
[X,Y]=np.meshgrid(longc,latgc)
for s in range (0,59): 
    if S.NAME[s] == str(STATE):
     print(S.NAME[s])
     multipolygon =Sgeo[s]
     print(len(multipolygon))
#for ij in range(len(multipolygon)):
in1=inpolygon(89,91,X,Y,multipolygon)
print(np.max(in1 ))  #,np.min(multipolygon[ij])
in11=in1  
print(np.max(in11))  

#-------------------------------------------------------------------------% 

# get fractional grid cells values of coarse cells

#-------------------------------------------------------------------------%
for I in range(1,IIPAR-1):
   for J in range(1,JJPAR-1):
      C = 0; dJ = (J-1),J,(J+1); dI = (I-1),I,(I+1); n = 2; err = 1
      if in1[J,I] == 0 and np.any(in1[dJ,dI]) == 1:
        C=1
      if in1[J,I] == 1 and np.any(in1[dJ,dI]) == 0:
        C=1    
      if C==0:
          continue
      while err >0.05:
          ln      = np.linspace((longc[I]-1/3),(longc[I]+1/3),n)
          lt      = np.linspace((latgc[J]-1/4),(latgc[J]+1/4),n)
          [X,Y]   = np.meshgrid(ln,lt)
          for ij in range(len(multipolygon)):
            tmp1 = np.array((inpolygon(len(ln),len(lt),X,Y,multipolygon[ij]))) 
          tmp=np.sum(tmp1)  
         # print(tmp)
          err     = abs(in1[J,I] - tmp / (n*n))
          in1[J,I] = tmp / (n*n)
          n       = n + 2 
        
STATE_MASK = (in1)
print(np.max(STATE_MASK))
#print(np.sum(STATE_MASK))
#-------------------------------------------------------------------------%
import csv
f = open(out_dir+'CITY_05x06.csv', 'w');writer = csv.writer(f);writer.writerows(CITY_MASK);f.close()
f = open(out_dir+'STATE_05x06.csv', 'w');writer = csv.writer(f);writer.writerows(STATE_MASK);f.close()


#%% (3) Calculate cost-function values

#% PM 2.5
NUM = 0; DEN = 0
for i in range(len(lo)):
    for j in range(len(la)):
       x=((abs(lonpm-lo[i])));x=np.where(x==np.min(x))# PM calculation for teh city only
       y=((abs(latpm-la[j])));y=np.where(y==np.min(y))
       xp=((abs(lonpo-lo[i])));xp=np.where(xp==np.min(xp))# for teh population over the city only
       yp=((abs(latpo-la[j])));yp=np.where(yp==np.min(yp))
       if np.isnan(PM[x,y]) == False and POP[yp,xp] >0: 
               NUM = NUM + (PM[x,y]) * (Z[j,i]) *(POP[yp,xp]) # Z is teh mask and calculation using pop data and O3
               DEN = DEN + (Z[j,i])  * (POP[yp,xp])
CFPM = NUM / DEN

#% O3
NUM = 0; DEN = 0
for i in range(len(lo)):
    for j in range(len(la)):
        x=((abs(longc-lo[i])));x=np.where(x==np.min(x))
        y=((abs(latgc-la[j])));y=np.where(y==np.min(y))
        xp=((abs(lonpo-lo[i])));xp=np.where(xp==np.min(xp))
        yp=((abs(latpo-la[j])));yp=np.where(yp==np.min(yp)) 
        if np.isnan(O3[x,y]) == False and POP[yp,xp] >0: 
               NUM = NUM + (O3[x,y]) * (Z[j,i]) *(POP[yp,xp])
               DEN = DEN + (Z[j,i])  * (POP[yp,xp])
CFO3 = NUM / DEN
        

#% NO2
NUM = 0; DEN = 0
for i in range(len(lo)):
    for j in range(len(la)):
        x=((abs(lonno-lo[i])));x=np.where(x==np.min(x))
        y=((abs(latno-la[j])));y=np.where(y==np.min(y))
        xp=((abs(lonpo-lo[i])));xp=np.where(xp==np.min(xp))
        yp=((abs(latpo-la[j])));yp=np.where(yp==np.min(yp))      
        if np.isnan(NO[x,y]) == False and POP[yp,xp] >0:
            NUM = NUM + ((NO[x,y])) * ((Z[j,i])) *(POP[yp,xp])
            DEN = DEN + (Z[j,i])  * (POP[yp,xp])     
CFNO = NUM / DEN
CFNO=CFNO.squeeze()



#%% (4) Emission scaling for gctm.gdt files        
for d in range(dmax):
    date=datetime(2010,11,30)+timedelta(days=d+1)
    if date.month <10:
        date_month1=str(date.month).rjust(1 + len(str(date.month)), '0')
    if date.month >9:
        date_month1=date.month   
    if date.day <10:
        date_day1=str(date.day).rjust(1 + len(str(date.day)), '0')
    if date.day >9:
        date_day1=date.day 
    vocfile = emis_dir+'ctm.VOC.'+str(date.year)+str(date_month1)+str(date_day1)+'.nc'
    noxfile = emis_dir+'ctm.NOx.'+str(date.year)+str(date_month1)+str(date_day1)+'.nc'
   # % NOx scaling
    if os.path.isfile(noxfile):
       ems=np.transpose(Dataset(noxfile, mode='r').variables['NOX-AN-S__NOx'][:])
       ems1 = np.sum(ems,axis=2)
       ems = np.sum(ems,axis=2)* efactor[0] #* np.transpose(GCM2)
       ems=np.multiply(ems,np.transpose(GCM2))# shorten 
       PMEMS[:,:,d,24] = np.divide(PMGDT[:,:,d,24] ,ems,out=np.zeros_like(PMEMS[:,:,d,24]), where=(PMGDT[:,:,d,24]!=0))
       O3EMS[:,:,d,24] = np.divide(O3GDT[:,:,d,24] , ems,out=np.zeros_like(O3EMS[:,:,d,24]), where=(O3GDT[:,:,d,24]!=0))
       NOEMS[:,:,d,24] = np.divide(NOGDT[:,:,d,24] , ems, out=np.zeros_like(NOEMS[:,:,d,24]), where=(NOGDT[:,:,d,24]!=0))
     #%VOC scaling for O3
    for k in range(k_o3):
        inames11=str(inames1[k]); kname=(inames11[12:len(inames11)-2].replace('_an',''));s11='ANTHSRCE__'+kname 
        try:
            temp=Dataset(vocfile, mode='r').variables[s11][:]
            temp=np.transpose(temp)
        except:
            continue
      
        if 'CO' in kname: 
                  ems=np.multiply(temp * efactor[1] , np.transpose(GCM2))
        else:
                 ems=np.multiply(temp * efactor[2],  np.transpose(GCM2))
        cf = np.mean(ems[ems>0]) / (1E1)        
        tmp = np.divide(O3GDT[:,:,d,k] , ems); tmp[ems < cf] = 0; O3EMS[:,:,d,k] = tmp     
#% remove infinities and NaNs from emission division
PMEMS[ np.isinf(PMEMS)] = 0;PMEMS[ np.isnan(PMEMS)] = 0
O3EMS[ np.isinf(O3EMS)] = 0;O3EMS[ np.isnan(O3EMS)] = 0
NOEMS[ np.isinf(NOEMS)] = 0;NOEMS[ np.isnan(NOEMS)] = 0



#%% (5) Apply cost-function scaling

#% PM2.5
sf    = CFPM / np.sum(PMGDT) 
PMEMS = PMEMS * sf

#% O3
sf    = CFO3 / np.sum(O3GDT) 
O3EMS = O3EMS * sf

#% NO2
sf    = CFNO / np.sum(NOGDT) 
NOEMS = NOEMS * sf

#%% (6) Calculate contributions
#% Read in NEI emissions

files = glob.glob(nei_dir+'nei_emis*.nc') 
E     = np.zeros((91,89,396,k_o3))
for f in range(len(files)):
    pathname = files[f]
    species=[]
    for i in range(k_o3):
        inames11=str(inames1[i])
        
        if '_an' in inames11: # why -an only?
                iname=(inames11[12:len(inames11)-2])
                if 'an1' in iname :   
                            iname=(iname.replace('_an1',''))
                if 'an2' in iname :   
                            iname=(iname.replace('_an1',''))            
                if 'PO' in iname :   
                    iname=''
                if 'PI' in iname :   
                    iname=(iname.replace('PI_an',''))
                if 'an' in iname :   
                    iname=(iname.replace('_an',''))

                if 'SO21' in iname:
                        iname='SO2'
     
                if 'an2' in iname:
                  iname=str('SO22')
                if 'NOX' in iname :
                   iname='NOx'
                if 'ISOP' in iname:
                   iname=''
        else:
            iname=''
        species.append(iname)
        if iname != '':       
         E[:,:,:,i] = E[:,:,:,i] + np.transpose(Dataset(pathname, mode='r').variables[iname][:])

        
        
#Calculate contributions

species_f0 = ["" for x in range(62)]
species_f1 = ["" for x in range(62)]
species_f2 = ["" for x in range(62)]
dJ_PM=   np.zeros((91,89,396,34))  
dJ_O3=   np.zeros((91,89,396,62))  
dJ_NO=   np.zeros((91,89,396,34))    
for i in range( 62):
    for p in range (0,3):
        if p==0:
              for k in range(len(inames0)):
                if inames1[i] == inames0[k]:
                    species_f0[k]=str(species[i])
                    dJ_PM[:,:,:,k] = E[:,:,:,i] * PMEMS[:,:,:,k]                   
        if p==1:
          for k in range(len(inames1)):
            if inames1[i] == inames1[k]:
                species_f1[k]=str(species[i])
                dJ_O3[:,:,:,k] = E[:,:,:,i] * O3EMS[:,:,:,k]
        if p==2:
          for k in range(len(inames2)):
            if inames1[i] == inames2[k]:
                species_f2[k]=str(species[i])
                dJ_NO[:,:,:,k] = E[:,:,:,i] * NOEMS[:,:,:,k]       
              
             
#%% (7) Preprocess health data

#% Load population data, mortality rates, and GBD codes
CPOP = pd.read_csv(ctypopfile,skiprows=1, names=['ID',	'Year',	'City',	'Country',	'Latitude',	'Longitude',	'Population',	'PM',	'Pw_PM',	'PAF_PM',	'Cases_PM',	'Rates_PM',	'O3',	'Pw_O3',	'PAF_O3',	'Cases_O3',	'Rates_O3',	'NO2',	'Pw_NO2',	'PAF_NO2',	'Cases_NO2',	'Rates_NO2',	'Pop_ped',	'GBDRegion',	'GBDSuperRegion',	'SDGRegion',	'WHORegion',	'WHOIncomeRegion',	'C40 Region'])
CITYPOP = np.array(CPOP[['Population']]).squeeze()
POP_PED = np.array(CPOP[['Pop_ped']]).squeeze()
Cind = np.array(CPOP[['City']]).squeeze()
for kkk in range(len(Cind)):
    if Cind[kkk] == NAME:
     ind1=kkk
POP_PED=POP_PED[ind1]
CITYPOP=CITYPOP[ind1]



GBDPOP = pd.read_csv(gbdpopfile,skiprows=1, names=['location_id',	'location_name',	'sex_id',	'sex_name',	'age_group_id',	'age_group_name	','year_id',	'measure_id',	'measure_name',	'metric_id','metric_name',	'val',	'upper',	'lower'])
Sind = np.array(GBDPOP[['location_name']]).squeeze()
Ssind=np.array(GBDPOP[['sex_name']]).squeeze()
ind1=np.where(Sind== STATE );ind1=np.array(ind1).squeeze()
GBDPOP=GBDPOP.iloc[ind1,:];Ssind=Ssind[ind1]
ind1=np.where(Ssind== SEX);ind1=np.array(ind1).squeeze()
GBDPOP=GBDPOP.iloc[ind1,:]

MORT  = pd.read_csv(mortfile,skiprows=1) ;ASTH =   pd.read_csv(asthfile,skiprows=1);A    = pd.read_csv(acodefile,skiprows=1) ;L    = pd.read_csv(lcodefile,skiprows=1) ;C    = pd.read_csv(ccodefile,skiprows=1) 




#% Process population data
POP = np.zeros((amax,1))
LB  = np.zeros((amax,1))
UB  = np.zeros((amax,1))
#% loop through age groups and calculate pop in age brackets
for a in range(amax):  
   aind =  np.where((GBDPOP.iloc[:,5] >= np.array(ages[a]))  & (GBDPOP.iloc[:,5] < np.array(ages[a])+5))#and (GBDPOP.iloc[:,5] < np.array(ages[a])+5))
   POP[a]=np.sum(GBDPOP.iloc[aind[0],11])
   LB[a]=np.sum(GBDPOP.iloc[aind[0],13])/np.sum(GBDPOP.iloc[aind[0],11]) 
   UB[a]=np.sum(GBDPOP.iloc[aind[0],12])/np.sum(GBDPOP.iloc[aind[0],11]) 
   if a == amax-1:
       temp11 = np.array(GBDPOP.iloc[:,11]);POP[a] = temp11[len(temp11)-1]
       temp11 = np.array(GBDPOP.iloc[:,13]);         LB[a] = temp11[len(temp11)-1] / POP[a]     
       temp11 = np.array(GBDPOP.iloc[:,12]);       UB[a] = temp11[len(temp11)-1] / POP[a]       
            
#% get all ages population
TOT     = np.sum(GBDPOP.iloc[:,11])   
TOT_LB  = np.sum(GBDPOP.iloc[:,13])   
TOT_UB  = np.sum(GBDPOP.iloc[:,12])   
#% calculate bounds of PED pop
POP_PED_LB = POP_PED * (TOT_LB -np.sum(POP*LB)) / (TOT-np.sum(POP))
POP_PED_UB = POP_PED * (TOT_UB -np.sum(POP*UB)) / (TOT-np.sum(POP))
#% calculate final city pop
POP  = CITYPOP * (POP / TOT)
POP_LB = LB * POP
POP_UB = UB * POP          


#% Process mortality rates
for i in range(len(L.iloc[:,0])):
    if L.iloc[i,2] == STATE:
        lid=(L.iloc[i,1])
rate      = np.zeros((amax,mmax))    
rate_lb   = np.zeros((amax,mmax))    
rate_ub   = np.zeros((amax,mmax))  

rate_asth = ASTH.iloc[np.where(ASTH.iloc[:,1]==lid)[0],7]
rate_asth_lb = ASTH.iloc[np.where(lid == ASTH.iloc[:,1])[0],9]
rate_asth_ub = ASTH.iloc[np.where(lid == ASTH.iloc[:,1])[0],8]
#% loop through mortality rates
for m in range(mmax):
   mid   = mortid[m]
#% loop through ages
   for a in range(amax):
    aid = ageids[a]
    inds = np.where((MORT.iloc[:,3] == aid) & (MORT.iloc[:,1] == lid) & (MORT.iloc[:,4] == mid))
    rate[a,m] = MORT.iloc[inds[0],7]
    rate_lb[a,m] = MORT.iloc[inds[0],9]
    rate_ub[a,m] = MORT.iloc[inds[0],8]
    
#%% (8) Health impact contribution

#% PM2.5
#%-------------------------------------------------------------------------%
#% load city-independent health data

rr_lb = np.transpose(Dataset(hlth_dir +'rr_lb.nc', mode='r').variables['rr_lb'][:])# what do they contain?
rr_ub = np.transpose(Dataset(hlth_dir +'rr_ub.nc', mode='r').variables['rr_ub'][:])
exp = np.transpose(Dataset(hlth_dir +'exp.nc', mode='r').variables['rexp'][:])
rr = np.transpose(Dataset(hlth_dir +'rr.nc', mode='r').variables['rr'][:])

#% initialize outputs
dJPMD    = np.zeros((IIPAR,JJPAR,dmax,k_pm))
dJPMD_lb = np.zeros((IIPAR,JJPAR,dmax,k_pm))
dJPMD_ub = np.zeros((IIPAR,JJPAR,dmax,k_pm))
AVAIL=np.zeros(np.shape(dJ_PM)[3])

#assign the caculations for the health facor calculations.
       
for ij in range(np.shape(dJ_PM)[3]):
 AVAIL[ij]  = (np.sum(dJ_PM[:,:,:,ij]))
#% calculate health impacts
for k in range(k_pm):
    if AVAIL[k] == 0: 
        continue 
    #print(k,AVAIL[k])
    amax = len(POP);    #% number of age groups
    omax = np.shape(rate)[1] #% number of health outcomes
    #print(omax)
    numVars = len(np.shape(dJ_PM[:,:,:,0]))#en(np.shape(dJ)); 
    sz=np.zeros((numVars))
    for n in range(numVars):
         sz[n]=np.shape(dJ_PM)[n]
    if numVars==2:
      D    = np.zeros((int(sz[0]),int(sz[1])))
      D_lb    = np.zeros((int(sz[0]),int(sz[1])))
      D_ub    = np.zeros((int(sz[0]),int(sz[1])))
    if numVars==3:
      D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
      D_ub    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
      D_lb   = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
    if numVars==4:
      D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
      D_lb    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
      D_ub    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
    base_pm= np.zeros((amax,omax))      
    base_pm_ub= np.zeros((amax,omax)) 
    base_pm_lb= np.zeros((amax,omax))
    for a in range(amax):#14,15
     for o in range(omax): 
         #print(a,o,np.shape(D))
         RR = rr[a,o,:] ;RR_lb = rr_lb[a,o,:];RR_ub = rr_ub[a,o,:]
         EXP = exp[a,o,:] 
         pop = POP[a];pop_lb = POP_LB[a];pop_ub = POP_UB[a]
         y0 = rate[a,o] / 1E5;y0_lb = rate_lb[a,o] / 1E5;y0_ub = rate_ub[a,o] / 1E5
         ind_1=[];ind_2=[]
         
         for ik in range (len (EXP)):
             if CFPM >= EXP[ik]: #???
              ind_1.append(ik)  
             if CFPM < EXP[ik]:
                  ind_2.append(ik)      
         ind_1=ind_1[len(ind_1)-1]
         ind_2=ind_2[0]
         z0    = EXP[ind_1]; zf    = EXP[ind_2] 
         frac  = (CFPM - z0)/(zf-z0) 
         RR_1 = RR[ind_1]; RR_2 = RR[ind_2]
         base_d1=np.divide(np.multiply(np.multiply(pop,y0),(RR_1 - 1)),RR_1);base_d1_lb=np.divide(np.multiply(np.multiply(pop_lb,y0_lb),(RR_1 - 1)),RR_1);base_d1_ub=np.divide(np.multiply(np.multiply(pop_ub,y0_ub),(RR_1 - 1)),RR_1)
         base_d2 = np.divide(np.multiply(np.multiply(pop,y0),(RR_2-1)),RR_2);base_d2_lb = np.divide(np.multiply(np.multiply(pop_lb,y0_lb),(RR_2-1)),RR_2);base_d2_ub = np.divide(np.multiply(np.multiply(pop_ub,y0_ub),(RR_2-1)),RR_2)
         base_d = base_d1 + frac * (base_d2 - base_d1);base_d_lb = base_d1_lb + frac * (base_d2_lb - base_d1_lb);base_d_ub = base_d1_ub + frac * (base_d2_ub - base_d1_ub)  
         #print(a,o,base_d)
         base_pm[a,o] = base_d 
         base_pm_ub[a,o] = base_d_ub
         base_pm_lb[a,o] = base_d_lb
         
         frac  = (-1)*(dJ_PM[:,:,:,k] + z0-CFPM)/(zf-z0)  
         pert_d = base_d1 + frac * (base_d2-base_d1);pert_d_ub = base_d1_ub + frac * (base_d2_ub-base_d1_ub) ;pert_d_lb = base_d1_lb + frac * (base_d2_lb-base_d1_lb)      
         D = D + base_d-pert_d 
         D_lb = D_lb + base_d_lb-pert_d_lb
         D_ub = D_ub + base_d_ub-pert_d_ub               
    dJPMD[:,:,:,k]=D
    dJPMD_ub[:,:,:,k]=D_ub
    dJPMD_lb[:,:,:,k]=D_lb
               
                


#-------------------------------------------------------------------------%

#% O3
#%-------------------------------------------------------------------------%
#% initialize outputs
dJO3D    = np.zeros((IIPAR,JJPAR,dmax,k_o3))
dJO3D_lb = np.zeros((IIPAR,JJPAR,dmax,k_o3))
dJO3D_ub = np.zeros((IIPAR,JJPAR,dmax,k_o3))
AVAIL=np.zeros(np.shape(dJ_O3)[3])
     
beta     = np.log(1.06) / 10
beta_lb  = np.log(1.03) / 10
beta_ub  = np.log(1.10) / 10
tmrel    = 32.4
tmrel_lb = 35.7
tmrel_ub = 29.1

for ij in range(np.shape(dJ_O3)[3]):
 AVAIL[ij]  = (np.sum(dJ_O3[:,:,:,ij]))
#% calculate health impacts
for k in range(k_pm):
     if AVAIL[k] == 0:  
        continue
     amax = len(POP);    #% number of age groups
     omax = 1 #% number of health outcomes
     numVars = len(np.shape(dJ_O3[:,:,:,0]))#en(np.shape(dJ)); 
     sz=np.zeros((numVars))
     for n in range(numVars):
          sz[n]=np.shape(dJ_O3)[n]
     if numVars==2:
       D    = np.zeros((int(sz[0]),int(sz[1])))
       D_lb    = np.zeros((int(sz[0]),int(sz[1])))
       D_ub    = np.zeros((int(sz[0]),int(sz[1])))
     if numVars==3:
       D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
       D_ub    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
       D_lb   = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
     if numVars==4:
       D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
       D_lb    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
       D_ub    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
     base_o3= np.zeros((amax,omax))      
     base_o3_ub= np.zeros((amax,omax)) 
     base_o3_lb= np.zeros((amax,omax))  
     for a in range(amax):
        for o in range(omax):
            #print(a,o)
            #% get age population
            pop = POP[a];
            #% get mortality rate for outcome and age group
            rateo3=rate[:,1]
            y0 = rateo3[a] / 1E5
            
            #% get relative risk from log-linear exposure response equation
            if CFO3 > tmrel:
              RR = np.exp(beta * (CFO3 - tmrel)) 
            else: 
              RR = 1  
            base_d =  np.divide(np.multiply(np.multiply(pop,y0),(RR-1)),RR) 
            #print(a,o,base_d)
            base_o3[a,o] = base_d
            temp = dJ_O3[:,:,:,k] 
            if CFO3 > tmrel:
                   RR = np.exp(beta * ((-1)* (temp + tmrel-CFO3))); 
            else:
                   RR = 1 

            pert_d = np.divide(np.multiply(np.multiply(pop , y0) , ( RR - 1 )) , RR);              
                       
               #% assign pert d to output array
            D = D + (base_d-pert_d); 
     dJO3D[:,:,:,k]=D

     
#%-------------------------------------------------------------------------%

#% NO2
#%-------------------------------------------------------------------------%

#% initialize outputs
dJNOD    = np.zeros((IIPAR,JJPAR,dmax,k_no))
dJNOD_lb = np.zeros((IIPAR,JJPAR,dmax,k_no))
dJNOD_ub = np.zeros((IIPAR,JJPAR,dmax,k_no))
AVAIL=np.zeros(np.shape(dJ_NO)[3])
     
beta     = np.log(1.26) / 10
beta_lb  = np.log(1.10) / 10
beta_ub  = np.log(1.37) / 10
tmrel    = 2
tmrel_lb = 5
tmrel_ub = 0

for ij in range(np.shape(dJ_NO)[3]):
 AVAIL[ij]  = (np.sum(dJ_NO[:,:,:,ij]))
#% calculate health impacts
for k in range(k_no):
     if AVAIL[k] == 0:  
        continue
     amax = np.size(POP_PED);    #% number of age groups
     omax = np.size(rate_asth) #% number of health outcomes
     numVars = len(np.shape(dJ_NO[:,:,:,0]))#en(np.shape(dJ)); 
     sz=np.zeros((numVars))
     for n in range(numVars):
          sz[n]=np.shape(dJ_NO)[n]
         
     if numVars==2:
       D    = np.zeros((int(sz[0]),int(sz[1])))
       D_lb    = np.zeros((int(sz[0]),int(sz[1])))
       D_ub    = np.zeros((int(sz[0]),int(sz[1])))
     if numVars==3:
       D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
       D_ub    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
       D_lb   = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
     if numVars==4:
       D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
       D_lb    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
       D_ub    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))
     base_no= np.zeros((amax,omax))      
     base_no_ub= np.zeros((amax,omax)) 
     base_no_lb= np.zeros((amax,omax)) 
     for a in range(amax):
        for o in range(omax):
            
            #% get age population
            pop = POP_PED#[a];
            #% get mortality rate for outcome and age group
            y0 = np.array(rate_asth)/ 1E5
            
            #% get relative risk from log-linear exposure response equation
            if CFNO > tmrel:
              RR = np.exp(beta * (CFNO - tmrel)) 
            else: 
              RR = 1 
            #print(a,o,RR)  
            base_d =  pop*y0*(RR-1)/(RR)#np.divide(np.multiply(np.multiply(pop,y0),(RR-1)),RR) 
            #print(a,o,base_d)
            base_no[a,o] = base_d
            temp = dJ_NO[:,:,:,k] 
            if CFNO > tmrel:
                   RR = np.exp(beta * ((-1)* (temp + tmrel-CFNO))); 
            else:
                   RR = 1 

            pert_d = np.divide(np.multiply(np.multiply(pop , y0) , ( RR - 1 )) , RR);              
                       
               #% assign pert d to output array
            D = D + (base_d-pert_d); 
     dJNOD[:,:,:,k]=D



#%-------------------------------------------------------------------------%

#% Calculate health scaling factors
PMHF = np.divide(dJPMD , dJ_PM); PMHF_LB = np.divide(dJPMD_lb , dJ_PM); PMHF_UB = np.divide(dJPMD_ub , dJ_PM)
O3HF = np.divide(dJO3D , dJ_O3); O3HF_LB = np.divide(dJO3D_lb , dJ_O3); O3HF_UB = np.divide(dJO3D_ub , dJ_O3)
NOHF = np.divide(dJNOD , dJ_NO); NOHF_LB = np.divide(dJNOD_lb , dJ_NO); NOHF_UB = np.divide(dJNOD_ub , dJ_NO)
#% Remove NaN values
PMHF[np.isnan(PMHF)] = 0; PMHF_LB[np.isnan(PMHF_LB)] = 0; PMHF_UB[np.isnan(PMHF_UB)] = 0
O3HF[np.isnan(O3HF)] = 0; O3HF_LB[np.isnan(O3HF_LB)] = 0; O3HF_UB[np.isnan(O3HF_UB)] = 0
NOHF[np.isnan(NOHF)] = 0; NOHF_LB[np.isnan(NOHF_LB)] = 0; NOHF_UB[np.isnan(NOHF_UB)] = 0


#%% (9) Write netCDF outputs for sensitivities

# %-----------------------
# % setup CAL2035 info
# %-----------------------

# % calculate date
now = datetime.now()
dt_string = now.strftime("%d%m%Y_%H%M%S")
import pandas as pd
from netCDF4 import Dataset,num2date,date2num
import datetime
# ----------------------
yy = '2010';mm='12';
date_series_min = pd.date_range(yy+'-'+mm, periods=396, freq='D') # This is DatetimeIndex, perhaps ordinary datetime objects are preferred...
datesout = date_series_min.to_pydatetime() 
unout = 'days since 2010-12-01 00:00:00'
ks = [k_pm,k_o3,k_no]     
import netCDF4 as nc
for p in range(3):    
#     % create file
    fn =out_dir+INIT+'_sens_'+pollutant[p]+'_'+dt_string+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    # ncout.createDimension('time',None);
    # timevar = ncout.createVariable('time','float32',('time'));timevar.setncattr('units',unout);
    # timevar[:]= date2num(datesout,unout);

    # % define dimensions
    time = ds.createDimension('time', dmax)
    lat = ds.createDimension('lat', JJPAR)
    lon = ds.createDimension('lon', IIPAR)
    
   # time = ds.createVariable('time', np.float64, ('time',))
    timevar = ds.createVariable('time','float32',('time'));timevar.setncattr('units',unout);
    lat = ds.createVariable('lat', 'f4', ('lat',))
    lon = ds.createVariable('lon', 'f4', ('lon',))
    
   # time[:]=pd.date_range(start="2010-12-01",end="2011-12-31")
    timevar[:]= date2num(datesout,unout);
    lat[:]=latgc
    lon[:]=longc
    for k in range(ks[p]):
        if p==0:
         spec=species_f0[k]
         if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(PMEMS[:,:,:,k])
        if p==1:
         spec=species_f1[k]  
         if spec != '':
             #print(p,spec)
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(O3EMS[:,:,:,k])
        if p==2:
            spec=species_f2[k]
            if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(NOEMS[:,:,:,k])
                
    ds.close()    
    
# %% (10) Write netCDF outputs for health fractions

# %-----------------------
# % setup CAL2035 info
# %-----------------------

for p in range(3):    
#     % create file
    fn =out_dir+INIT+'_health_'+pollutant[p]+'_'+dt_string+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    
    # % define dimensions
    time = ds.createDimension('time', dmax)
    lat = ds.createDimension('lat', JJPAR)
    lon = ds.createDimension('lon', IIPAR)
    
    timevar = ds.createVariable('time','float32',('time'));timevar.setncattr('units',unout);
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))

    timevar[:]= date2num(datesout,unout);
    lats[:]=latgc
    lons[:]=longc
    for k in range(ks[p]):
        if p==0:
         spec=species_f0[k]
         if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(PMHF[:,:,:,k])
        if p==1:
         spec=species_f1[k]  
         if spec != '':
            # print(p,spec)
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(O3HF[:,:,:,k])
        if p==2:
            spec=species_f2[k]
            if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(NOHF[:,:,:,k])
                
    ds.close()  
    
    
# %% (11) Write netCDF outputs for health fractions LOwer bounds 

# %-----------------------
# % setup CAL2035 info
# %-----------------------

for p in range(3):    
#     % create file
    fn =out_dir+INIT+'_health_lb_'+pollutant[p]+'_'+dt_string+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    
    # % define dimensions
    time = ds.createDimension('time', dmax)
    lat = ds.createDimension('lat', JJPAR)
    lon = ds.createDimension('lon', IIPAR)
    
    timevar = ds.createVariable('time','float32',('time'));timevar.setncattr('units',unout);
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))

    timevar[:]= date2num(datesout,unout);
    lats[:]=latgc
    lons[:]=longc
    for k in range(ks[p]):
        if p==0:
         spec=species_f0[k]
         if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(PMHF_LB[:,:,:,k])
        if p==1:
         spec=species_f1[k]  
         if spec != '':
            # print(p,spec)
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(O3HF_LB[:,:,:,k])
        if p==2:
            spec=species_f2[k]
            if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(NOHF_LB[:,:,:,k])
                
    ds.close() 


# %% (12) Write netCDF outputs for health fractions UPPER bounds 

# %-----------------------
# % setup CAL2035 info
# %-----------------------

for p in range(3):    
#     % create file
    fn =out_dir+INIT+'_health_ub_'+pollutant[p]+'_'+dt_string+'.nc'
    ds = nc.Dataset(fn, 'w', format='NETCDF4')
    
    # % define dimensions
    time = ds.createDimension('time', dmax)
    lat = ds.createDimension('lat', JJPAR)
    lon = ds.createDimension('lon', IIPAR)
    
    timevar = ds.createVariable('time','float32',('time'));timevar.setncattr('units',unout);
    lats = ds.createVariable('lat', 'f4', ('lat',))
    lons = ds.createVariable('lon', 'f4', ('lon',))

    timevar[:]= date2num(datesout,unout);
    lats[:]=latgc
    lons[:]=longc
    for k in range(ks[p]):
        if p==0:
         spec=species_f0[k]
         if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(PMHF_UB[:,:,:,k])
        if p==1:
         spec=species_f1[k]  
         if spec != '':
             print(p,spec)
             value = ds.createVariable(str(spec), 'f4', ('time','lat','lon'))
             value[:,:,:]=np.transpose(O3HF_UB[:,:,:,k])
        if p==2:
            spec=species_f2[k]
            if spec != '':
             value = ds.createVariable(str(spec), 'f4', ( 'time','lat','lon'))
             value[:,:,:]=np.transpose(NOHF_UB[:,:,:,k])
                
    ds.close()          
    
    
#print(date(2014, 7, 2) )   
    
# leave for now for teh file info in the nc file
#     % create file
# ncid = netcdf.create(fname_out,'netcdf4') 

# % define dimensions
# latid    = netcdf.defDim(ncid,'lat',JJPAR)
# lonid    = netcdf.defDim(ncid,'lon',IIPAR)  
# dayid    = netcdf.defDim(ncid,'time',dmax)

# % define 1d variables
# lonvarid = netcdf.defVar(ncid,'lon','NC_DOUBLE',lonid)
# netcdf.putAtt(ncid,lonvarid,'units','degrees_east')
# netcdf.putAtt(ncid,lonvarid,'long_name','longitude')
# netcdf.putAtt(ncid,lonvarid,'comment','centre of grid cell')

# latvarid = netcdf.defVar(ncid,'lat','NC_DOUBLE',latid)
# netcdf.putAtt(ncid,latvarid,'units','degrees_north')
# netcdf.putAtt(ncid,latvarid,'long_name','latitude')
# netcdf.putAtt(ncid,latvarid,'comment','centre of grid cell')   

# dayvarid = netcdf.defVar(ncid,'time','NC_DOUBLE',dayid)
# netcdf.putAtt(ncid,dayvarid,'units','days since 1981-01-01 00:00:00')
# netcdf.putAtt(ncid,dayvarid,'long_name','reference time of sst field')
# netcdf.putAtt(ncid,dayvarid,'comment','')  
# netcdf.putAtt(ncid,dayvarid,'axis','T')







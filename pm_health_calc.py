#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:37:06 2022

@author: schakrab
"""

def pm_health_calc(pop,rate1,exps1,rr1,cfn,dJ):




 import numpy as np
 amax = len(pop);    #% number of age groups
 omax = np.shape(rate1)[1] #% number of health outcomes

 numVars = len(np.shape(dJ))#en(np.shape(dJ)); 
 sz=np.zeros((numVars))

 for n in range(numVars):
     sz[n]=np.shape(dJ)[n]
 
#% Initializations
 if numVars==2:
  D    = np.zeros((int(sz[0]),int(sz[1])))
 if numVars==3:
     D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])))
 if numVars==4:
      D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])))    

 base = np.zeros((amax,omax));

 for a in range(14,15):#14,15
    for o in range(5,6): 
        print(a,o)
        RR = rr1[a,o,:] 
        EXP = exps1[a,o,:] 
        POP = pop[a]
        y0 = rate1[a,o] / 1E5
        ind_1=[];ind_2=[]
        for ik in range (len (EXP)):
            if cfn >= EXP[ik]:
             ind_1.append(ik)  
            if cfn < EXP[ik]:
                 ind_2.append(ik)      
        ind_1=ind_1[len(ind_1)-1]
        ind_2=ind_2[0]
        z0    = EXP[ind_1]; zf    = EXP[ind_2] 
        frac  = (cfn - z0)/(zf-z0) 

        RR_1 = RR[ind_1]; RR_2 = RR[ind_2]
        base_d1=np.divide(np.multiply(np.multiply(POP , y0) , ( RR_1 - 1 )) , RR_1)
        base_d2 = np.divide(np.multiply(np.multiply(POP , y0) , ( RR_2 - 1 )) , RR_2)

        base_d = base_d1 + frac * (base_d2 - base_d1)     

        base[a,o] = base_d;

        if numVars == 2:
         for i in range(int(sz[0])):
            for j in range(int(sz[1])):
               temp = dJ[i,j];
               frac  = (cfn - temp - z0)/(zf-z0) 
               pert_d = base_d1 + frac * (base_d2-base_d1)      
               D[i,j] = D[i,j] + base_d-pert_d                
          
        elif numVars == 3:
         for i in range(int(sz[0])):
            for j in range(int(sz[1])):
              for kk in range(int(sz[2])):
               frac  = (cfn - dJ[i,j,kk] - z0)/(zf-z0)    
               pert_d = base_d1 + frac * (base_d2-base_d1)
               D[i,j,kk] = D[i,j,kk] + base_d-pert_d      
        elif numVars == 4:
         for i in range(int(sz[0])):
          for j in range(int(sz[1])):
            for kk in range(int(sz[2])):
             for l in range(int(sz[3])):  
               temp = dJ[i,j,kk,l]
               frac  = (cfn - temp - z0)/(zf-z0)    
               pert_d = base_d1 + frac * (base_d2-base_d1)      
               D[i,j,kk,l] = D[i,j,kk,l] + [base_d-pert_d]        
 #print(D)       
 return(base)
 return(D)     
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 09:49:59 2022

@author: schakrab
"""

def no_health_calc(pop,rate1,cfn,dJ,beta,tmrel):
# % Function "o3_health_calc" takes inputs of population data (pop),
# % mortality rates (rate), a base cost-function value (cfn), and an adjoint 
# % contribution array (dJ) and calculates the total and source apportionment of
# % heatlh impacts.

# % Variables:
# %-------------------------------------------------------------------------%
# % pop  - m x 1 vector of age-stratified population data
# % rate - m x n array of outcome and age seperated mortality rates
# % cfn  - scalar cost-function value of o3exposure
# % dJ   - 2D, 3D, or, 4D array of sources that contribute to the cost-function 
# %        exposure

# % Outputs
# %-------------------------------------------------------------------------%
# % base - total o3 exposure deaths from cost-function
# % D    - an array of source apportionment deaths

# %% Set-up

# % Assignments
 import numpy as np
 amax = len(pop);    #% number of age groups
 omax = np.shape(rate1)[1]#;; % number of health outcomes
 #% get size of source apportionment array
 numVars = len(np.shape(dJ))#en(np.shape(dJ)); 
 sz=np.zeros((numVars))
 for n in range(numVars):
     sz[n]=np.shape(dJ)[n]
#% Initializations
 if numVars==2:
  D    = np.zeros((int(sz[0]),int(sz[1])));
 if numVars==3:
     D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2])));
 if numVars==4:
      D    = np.zeros((int(sz[0]),int(sz[1]),int(sz[2]),int(sz[3])));    
 base = np.zeros((amax,omax));

#%% Calculate premature deaths

#% loop through age groups and outcomes
 for a in range(amax):
    for o in range(omax):
         
        #% get age population
        POP = pop[a];
        #% get mortality rate for outcome and age group
        y0 = rate1[a,o] / 1E5;
        #% get relative risk from log-linear exposure response equation
        if cfn >= tmrel:
          RR = np.exp(beta * (cfn - tmrel)); 
        else: 
          RR = 1
                        
        base_d = np.divide(np.multiply(np.multiply(POP , y0) , ( RR - 1 )) , RR);             
        
        #% assign base deaths to output array
        base[a,o] = base_d;

        
        #% loop through source groups depending on number of variables
        if numVars == 2:
         for i in range(int(sz[0])):
            for j in range(int(sz[1])):
               #% get source contribution
               temp = dJ[i,j]; 
               
               #% get relative risk from log-linear exposure response equation
               if cfn > tmrel:
                   RR = np.exp(beta * (cfn - temp - tmrel)); 
               else:
                   RR = 1; 
                
               
               #% calculate cost-function premature deaths for both risks
               pert_d = np.divide(np.multiply(np.multiply(POP , y0) , ( RR - 1 )) , RR);            
                       
               #% assign pert d to output array
               D[i,j] = D[i,j] + (base_d-pert_d);                
           
        if numVars == 3:
          for i in range(int(sz[0])):
           for j in range(int(sz[1])):
            for k in range(int(sz[2])):
               #% get source contribution
               temp = dJ[i,j,k]; 
               
               #% get relative risk from log-linear exposure response equation
               if cfn > tmrel:
                   RR = np.exp(beta * (cfn - temp - tmrel)); 
               else:
                   RR = 1 
                        
               
               #% calculate cost-function premature deaths for both risks
               pert_d = np.divide(np.multiply(np.multiply(POP , y0) , ( RR - 1 )) , RR);   
               
               #% assign pert d to output array
               D[i,j,k] = D[i,j,k] + (base_d-pert_d);                

        if numVars == 4:
            for i in range(int(sz[0])):
                for j in range(int(sz[1])):
                    for k in range(int(sz[2])):
                        for l in range(int(sz[3])):  
               #% get source contribution
                           temp = dJ[i,j,k,l]; 
               
               #% get relative risk from log-linear exposure response equation
                           if cfn > tmrel:
                                   RR = np.exp(beta * (cfn - temp - tmrel)); 
                           else:
                                   RR = 1; 
      
               
               #% calculate cost-function premature deaths for both risks
                           pert_d = np.divide(np.multiply(np.multiply(POP , y0) , ( RR - 1 )) , RR);    
               
               #% assign pert d to output array
                           D[i,j,k,l] = [(i,j,k,l)]+ (base_d-pert_d);                

 return(base)
 return(D) 

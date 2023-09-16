#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[]:

"""    

Script to find melting point using fsolve.
Functional form of free energies using Hessian and Equilibrium composition 
data from pandat
free energies expressed as second order to compositions

"""

# In[]:

import pandas as pnd
from pandas import DataFrame as df
import numpy as np
from scipy.interpolate import UnivariateSpline
import re
from scipy.optimize import fsolve 

# In[]:

"""    

Inputs for the script are taken from text or csv files. 
Compositions and hessians for each phase are taken as input  

"""

comp_filename  = "Composition_FCC_A1.csv"    

with open(comp_filename, 'r') as my_file:
    
    flines=my_file.readlines()

    for line in flines:
        searchstatus = re.match(r'\bTemp\b', line)
      
        if searchstatus:
            
            line = line.split(",")
            num_components = int(0.5*(len(line) - 1)) # num_components : Independent Ones
            break

comp_df = df(pnd.read_csv(comp_filename,delimiter=","))
comp_data = comp_df.sort_values(by="Temp")

T_array      = comp_data[comp_data.columns[0]].to_numpy()

# In[]:
 
# order of composition : Solid for all Independent components and then liquid   
# composition of Liquid and Solid are written as functions of temperature
# The composition of 1st element can be found from function Cs_fun[1](T)    
# and Cl_fun[1](T)   
 
Cs_fun = {}
Cl_fun = {}    

for i in range(num_components):
    
    Sol_Col_Num = i+1
    Cs_array = comp_data[comp_data.columns[Sol_Col_Num]].to_numpy()
    Cs_fun[Sol_Col_Num] = UnivariateSpline(T_array, Cs_array, k=3, s=0)

Col_Num = Sol_Col_Num
    
for i in range(num_components):
    
    Liq_Col_Num = i+1  
    Cl_array = comp_data[comp_data.columns[Col_Num+1+i]].to_numpy()
    Cl_fun[Liq_Col_Num] = UnivariateSpline(T_array, Cl_array, k=3, s=0)  

# In[]:
    
#print(Cs_fun[1](1560))    

# In[]:

# HSN is the hessian from csv file. Functions for ternary As_fun[11](T), 
# As_fun[22](T), As_fun[33](T), As_fun[12](T),As_fun[13](T), As_fun[23](T)     

HSN_solid_filename = "HSN_FCC_A1.csv"
 
HSN_sol_df = df(pnd.read_csv(HSN_solid_filename,delimiter=","))
HSN_sol_data = HSN_sol_df.sort_values(by="Temp")

As_fun = {}

for i in range(num_components):

    arr_Nam = str(i+1) + str(i+1)
    arr_Num = int(arr_Nam)
    As_arr = 0.5*HSN_sol_data[HSN_sol_data.columns[i+1]].to_numpy()
    As_fun[arr_Num] = UnivariateSpline(T_array, As_arr, k=3, s=0)
    
Col_Num = num_components
    
for i in range(1,num_components+1,1):
    for j in range(i+1,num_components+1,1):
    
        Col_Num += 1
        arr_Nam = str(i) + str(j)
        arr_Num = int(arr_Nam)
        As_arr = HSN_sol_data[HSN_sol_data.columns[Col_Num]].to_numpy()
        As_fun[arr_Num] = UnivariateSpline(T_array, As_arr, k=3, s=0) 

# In[]:

#print(As_fun[12](1560))    
            
# In[]:

# HSN is the hessian from csv file. Functions for ternary Al_fun[11](T), 
# Al_fun[22](T), Al_fun[33](T), Al_fun[12](T),Al_fun[13](T), Al_fun[23](T)     
  
HSN_liq_filename = "HSN_LIQUID.csv"
 
HSN_liq_df = df(pnd.read_csv(HSN_liq_filename,delimiter=","))
HSN_liq_data = HSN_liq_df.sort_values(by="Temp")

Al_fun = {}

for i in range(num_components):

    arr_Nam = str(i+1) + str(i+1)
    arr_Num = int(arr_Nam)
    Al_arr = 0.5*HSN_liq_data[HSN_liq_data.columns[i+1]].to_numpy()
    Al_fun[arr_Num] = UnivariateSpline(T_array, Al_arr, k=3, s=0)
    
Col_Num = num_components
    
for i in range(1,num_components+1,1):
    for j in range(i+1,num_components+1,1):
    
        Col_Num += 1
        arr_Nam = str(i) + str(j)
        arr_Num = int(arr_Nam)
        Al_arr = HSN_liq_data[HSN_liq_data.columns[Col_Num]].to_numpy()
        Al_fun[arr_Num] = UnivariateSpline(T_array, Al_arr, k=3, s=0)

# In[]:

#print(Al_fun[12](1560))    

# In[]:

# B is calculated for only for solid phase for each independent component from   
# hessians and equilibrium calculations. B_fun[1](T), B_fun[2](T) .....

B_fun = {}    

for j in range(1,num_components+1,1):
    
    B_Nm = str(j)
    B_Nu = int(B_Nm)
    B_arr = np.zeros(len(T_array))
    
    for n in range(len(T_array)):     
        T = T_array[n]
        
        for i in range(1,num_components+1,1):
           
            A_Nm = str(i) + str(j)
            A_Num = int(A_Nm)
            
            if i == j:
            
                B_arr[n] += 2*(Al_fun[A_Num](T)*Cl_fun[B_Nu](T) \
                     - As_fun[A_Num](T)*Cs_fun[B_Nu](T)) 
                    
            if i != j:
               
                if i > j:
                   
                   A_Num = int(str(j)+ str(i))

                B_arr[n] +=  Al_fun[A_Num](T)*Cl_fun[int(str(i))](T) \
                     - As_fun[A_Num](T)*Cs_fun[int(str(i))](T)    
              
    B_fun[B_Nu] = UnivariateSpline(T_array, B_arr, k=3, s=0)

# In[]:

#print(B_fun[2](1655))
 
# In[]:

# C is calculated for only for solid phase for each independent component from   
# hessians and equilibrium calculations. C_fun(T)
    
C_arr = np.zeros(len(T_array)) 
   
for n in range(len(T_array)):
    
    T = T_array[n]
    
    for i in range(1,num_components+1,1):
        for j in range(i,num_components+1,1): 
             
            A_Nm = str(i) + str(j)
            A_Num = int(A_Nm)
            
            C_arr[n] += As_fun[A_Num](T)*Cs_fun[int(str(i))](T)*Cs_fun[int(str(j))](T) \
                        - Al_fun[A_Num](T)*Cl_fun[int(str(i))](T)*Cl_fun[int(str(j))](T)

    
C_fun = UnivariateSpline(T_array, C_arr, k=3, s=0)   
 
# In[]:

#print(C_fun(1625))    

# In[]:

def fun_solid(T,c):    

    """    
     c is of Size num_component(Independent) and has components
     T is Temperature
     
     returns free energy for Solid Phase 
    """
    
    f_sol = 0
    
    f_sol += C_fun(T)
    
    for i in range(1, num_components + 1,1):
        
        for j in range(i, num_components +1,1):
            
            A_Nm = str(i) + str(j)
            A_Num = int(A_Nm)
            
            f_sol += As_fun[A_Num](T)*c[i-1]*c[j-1]   
        
        f_sol += B_fun[int(str(i))](T)*c[i-1] 
   
    return f_sol
# In[]:
    
def fun_liq(T,c):    
    
    """    
     c is of Size num_component(Independent) and has components
     T is Temperature
     
     returns free energy for Liquid Phase 
    """
    
    f_liq = 0
    
    for i in range(1, num_components + 1,1):
        for j in range(i, num_components +1,1):
           
            A_Nm = str(i) + str(j)
            A_Num = int(A_Nm)
             
            f_liq += Al_fun[A_Num](T)*c[i-1]*c[j-1]              
     
    return f_liq  

# In[]:
    
def equations_solve(c_T_arr,c_l):    
 
    """
    equations to solve for Melting Temperature and Solid_Composition 
    for all Components when Liquid Composition is known
                         
    """
    
    T = c_T_arr[-1]
    c_s = c_T_arr[0:-1]
    
   # print(T,c_s)
    
    eq_ua = []    
     
    eq_1 = fun_liq(T,c_l) - fun_solid(T,c_s)
    
    for i in range(1,num_components+1,1):
        
        sum_s = 0
        sum_l = 0
        
        sam_l = 2*Al_fun[int(str(i)+ str(i))](T)*c_l[i-1] 
        sam_s = 2*As_fun[int(str(i)+ str(i))](T)*c_s[i-1]
        
        sam_s = sam_s + B_fun[int(str(i))](T)
        
        for j in range(1,num_components+1,1): 
            
            if j != i:
               
               fun_nam = int(str(i)+ str(j)) 
               
               if j < i:               
                   
                   fun_nam = int(str(j)+ str(i))
                   
               sum_l += Al_fun[fun_nam](T)*c_l[j-1]
               sum_s += As_fun[fun_nam](T)*c_s[j-1]
        
        sum_l += sam_l
        sum_s += sam_s
        eq_ua.append(sum_l - sum_s)
        eq_1 += sum_s*c_s[i-1] - sum_l*c_l[i-1]
        
    eq_ua.append(eq_1)
    eq_ua = np.array(eq_ua)      
  
    return eq_ua

# In[]:

def find_equilibria(Guess,c_l):
    
    """
    Find melting point of liquid with composition c_l and composition of solid
    in equilibrium with c_l
    
    """
    
    res = fsolve(equations_solve,Guess, args = c_l)
    
    err = equations_solve(res, c_l)
    
    if np.max(err) > 1e-4:
    
         print(' error :',  err )
    
    return  res

# In[]:

T_melt =  1590

Cl_melt = np.zeros(num_components) 
Cs_melt = np.zeros(num_components) 

for i in range(num_components):
    
    Cl_melt[i] = Cl_fun[int(str(i+1))](T_melt) 
    Cs_melt[i] = Cs_fun[int(str(i+1))](T_melt)

# In[]:

Guess_melt = np.append(Cs_melt, T_melt)    

c_liq =  [0.101734,0.134514]   

res_ult = find_equilibria(Guess_melt, c_liq)
c_sol = res_ult[0:-1]
T_found = res_ult[-1]

print( 'Liq composition : ', c_liq )
print( 'Sol composition : ', c_sol )
print( 'Equilibrium Temp  : ', T_found)



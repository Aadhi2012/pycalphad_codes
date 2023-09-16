#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[]:
    
import numpy as np
from pycalphad import Database, equilibrium
import pycalphad.variables as v
import Liquidus_Temp as Solv_Liq_Temp 
from operator import itemgetter
import matplotlib.pyplot as plt 
import VR_tip as VR_TIP
from scipy.interpolate import interp1d
from pandas import DataFrame as df

# In[]:

db_sys = Database("alzn_mey.tdb")

phases_sys = list(db_sys.phases.keys())
phases_sys = sorted(phases_sys) 
print('phase Present in system are', phases_sys)
print('')

# In[]:

#Solid_Phase = sys.argv[2]  # Input no. 2
#Liquid_Phase = sys.argv[3] # Input no. 3

Solid_Phase = "FCC_A1"
Liquid_Phase = "LIQUID"

phases_calc = [Solid_Phase,Liquid_Phase]  # Phases for Calculation  
phases_calc_sort = sorted(phases_calc)    # Phases Sorted       

print('Solid Phase :', phases_calc[0])
print('Liquid Phase :', phases_calc[1])
print('')

# In[]:

component_sys = list(db_sys.elements)
for component in component_sys:

    if component == '/-':
        component_sys.remove(component)               

print('Components in system are', component_sys)
print('')

# In[]:
    
calc_component_sys = [x for x in component_sys if x != 'VA']
calc_component_sys.sort()

#Compo_calc_1 = sys.argv[4]  # Component for Calculation  Input no. 4

Compo_calc_1 = calc_component_sys[0]  # Component for Calculation
ind_ex_1 = calc_component_sys.index(Compo_calc_1)    

print('Component for calculation :' , Compo_calc_1)
print('')

# In[]:

def equi_composition_phase(Temperature,Composition,phase,phases_in_sys):  
    
    
    """
     Function to compute equilibrium phase composition at a Temperature 
 
     Inputs:
 
     Temperature
     Composition : Alloy Composition    
     phase : phase for composition is found (ex: "LIQUID")  
     phases_in_sys : list of phases for the calculation (ex: ["FCC_A1","LIQUID"])  
 
     return:
         
     Output array containing equilibrium composition for phase    
     
    """
    phase = [phase]
    
    #print(phase)
    
    eq_comp = equilibrium(db_sys, component_sys, phases_in_sys,{v.X(Compo_calc_1):Composition,v.T:Temperature,v.P: 101325},output='GM')
    Phase_at_equi = np.array(eq_comp.Phase.squeeze())
    Phase_at_equi = [x for x in Phase_at_equi if x] # Eliminate Nan Entries
    Phase_at_equi = sorted(Phase_at_equi)
     
    # print('Phase_at_equilibrium = ',Phase_at_equi)
     
    if phase[0] in Phase_at_equi:
         
        phase_compo = np.array(eq_comp.X.where(eq_comp.Phase == phase).sel(component = Compo_calc_1).squeeze())
        phase_compo = phase_compo[~(np.isnan(phase_compo))]
        return phase_compo 
     
    else:
         
        print('Phase not present for the Temperature and Composition')
        return
    
# In[]:
    
#print(equi_composition_phase(600, 0.4, Solid_Phase, phases_sys))
# In[ ]:

def phase_in_equilibrium(Temperature,Composition,phases_in_sys):
    
    """
    
    Function to computes equilibrium phases(sorted) at a Temperature 
    for a given composition 
    
    Input :
    
    Temperature
    Composition :Alloy Composition     
    phases_in_sys:list of phases for the calculation (ex: ["FCC_A1","LIQUID"])  
 
    return :
         
    Output array containing equilibrium phases(sorted)    
    
    """
    eq_com = equilibrium(db_sys, component_sys,phases_in_sys,{v.X(Compo_calc_1):Composition,v.T:Temperature, v.P: 101325},output='GM')
    
    Phase_at_equi = np.array(eq_com.Phase.squeeze())
    Phase_at_equi = [x for x in Phase_at_equi if x]
    Phase_at_equi = sorted(Phase_at_equi)
    
  #  print('Phases at equilibrium for Temperature Composition :', Phase_at_equi)
    
    return Phase_at_equi

# In[]:
    
#print(phase_in_equilibrium(600, 0.4, phases_sys)) 


# In[]:
    
c_alloy = 0.8168 

print(' Alloy Composition :',c_alloy)

T_resl = Solv_Liq_Temp.Temp_Range_Equilibrium(c_alloy, phases_sys, 'FCC_A1')
Melt_Temp = T_resl[1]
Solidus_Temp = T_resl[0]

Tm_C0 = Melt_Temp 

slope = Solv_Liq_Temp.sign_slope_Liquidus(c_alloy,Tm_C0,phases_calc_sort)

# In[]:

if slope > 0:

    c_perturb = 0.0005

else:

    c_perturb = -0.0005
print(c_perturb)

# In[]:
    
del_c = 0.01    

Compo_LiqT_arr = VR_TIP.Liquidus_Tem_entire(c_alloy,Tm_C0,del_c)    
Compo_LiqT_arr = np.array(Compo_LiqT_arr)
Compo_LiqT_arr = sorted(Compo_LiqT_arr, key=itemgetter(0))
Compo_LiqT_arr = np.array(Compo_LiqT_arr)

# In[]:

Compo_Liq_arr = Compo_LiqT_arr[:,0]
T_Liq = Compo_LiqT_arr[:,1]

# In[]:
Compo_Sol_arr = np.zeros(len(Compo_Liq_arr))
 
for i in range(len(Compo_Liq_arr)):
 
    Compo_Sol_arr[i] = equi_composition_phase(T_Liq[i],Compo_Liq_arr[i]+c_perturb,Solid_Phase,phases_sys)[0]

# In[]:

f_cubic_Liquidus = interp1d(Compo_Liq_arr, T_Liq, kind='cubic') # Liquidus Temp as a function of Liq Compo
f_cubic_Liq_Sol = interp1d(Compo_Liq_arr, Compo_Sol_arr, kind='cubic') 

f_cubic_Sol_Liq = interp1d(Compo_Sol_arr, Compo_Liq_arr, kind='cubic') 
f_cubic_Solid_Compo = interp1d(T_Liq,Compo_Sol_arr, kind='cubic') # Solid Compo as a function of Liquidus Temp
f_cubic_Liq_Compo = interp1d(T_Liq,Compo_Liq_arr, kind='cubic') # Solid Compo as a function of Liquidus Temp

# In[]:

df_s = 0.001    

# In[]:

def new_liq(c_liq,c_liq_eq,c_sol_eq,f_s):

    d_cl = (c_liq_eq - c_sol_eq)*(df_s)/(1-f_s)
    
    cl_new = c_liq + d_cl
    
    return cl_new

# In[]:    

def scheil_solid(c_alloy):

    T_start = Tm_C0
    
    f_s = 0
    c_liq = f_cubic_Liq_Compo(Tm_C0).item()
    f_l_arr , c_liq_arr,T_Liq_arr_sche = [], [], []
    T_current = T_start
    
    c_liq_arr.append(c_liq)
    f_l_arr.append(1 - f_s)
    T_Liq_arr_sche.append(T_current)

    while f_s < 0.6:
        
          c_liq_eq = f_cubic_Liq_Compo(T_current).item()
          c_sol_eq = f_cubic_Solid_Compo(T_current).item()
          c_liq_new = new_liq(c_liq, c_liq_eq,c_sol_eq ,f_s)             
       #   print(c_liq,c_liq_new)
    #      phases_current = phase_in_equilibrium(T_current, c_liq + c_perturb, phases_sys)
          c_liq = c_liq_new
          T_current = f_cubic_Liquidus(c_liq).item()
          f_s += df_s
          c_liq_arr.append(c_liq)
          f_l_arr.append(1-f_s)
          T_Liq_arr_sche.append(T_current)  
        #  print(f_s)
          
    return c_liq_arr, f_l_arr, T_Liq_arr_sche 
    
# In[]:
    
arr_A, arr_B, arr_c = scheil_solid(c_alloy)    

df_scheil = df({'c_liq': arr_A, 'f_l': arr_B, 'T': arr_c})
df_scheil.to_csv('scheil',index = False)

# In[]:
    
plt.plot(arr_B,arr_A)    


#!/usr/bin/env python
# coding: utf-8

# In[1]:

import numpy as np
import sys
from pycalphad import Database, equilibrium
import pycalphad.variables as v
import scipy.special as sp
import scipy.optimize as opti
from scipy.interpolate import interp1d
import re
import pandas as pnd
from pandas import DataFrame as df

import Liquidus_Tem as Solv_Liq_Temp
import Gibb_Thomson_Coef as GTC
import Liquidus_Tem_Entire as Liq_Tem_ent

# In[47]:

"""  

 This script is used to find Velocity and Radius of dendritic tip using LGK
 formulation for a given undercooling.
 
 Functions:
 
 find_V_R_tip() : finds Velocity and radius of dendritic tip  
 VR_entire_range() : finds Velocity and radius of dendritic tip for entire 
                     solidification range  

"""    
# In[]: 
"""    

Inputs read from .in input file taken as argument using sys.argv

Inputs read from input file:
    
tdbfname, NUMCOMPONENTS, COMPONENTS, tdb_phases, num_thermo_phases, ceq    

"""
print('Reading Input file')

#with open(sys.argv[1], 'r') as my_file:    
with open('Input.in', 'r') as my_file:    
    
    flines=my_file.readlines()

    flg_tdbf, flg_Comps, flg_Num_Comps, flg_Phase, flg_Num_Phases = 0,0,0,0,0  
    
    for line in flines:
       
        searchstatus = re.match(r'\btdbfname\b', line)
        if searchstatus:
            flg_tdbf = 1
            tdbfline = line
            tdbfname = tdbfline.split("=")[1].replace(";", "").strip()
            #print("Tdb file name : ", tdbfname)

        searchstatus = re.match(r'\bNUMCOMPONENTS\b', line)
        if searchstatus:
            flg_Num_Comps = 1
            N_COMPONENTSline = line
            compts = N_COMPONENTSline.split("=")[1].replace(';', '').strip()
            N_comps = int(np.array(re.sub(r"[{} ]","",compts).split(","))[0])
            #print("Num of Comps : ",N_comps)
        
        searchstatus = re.match(r'\bCOMPONENTS\b', line)
        if searchstatus:
            flg_Comps = 1
            COMPONENTSline = line
            components1 = COMPONENTSline.split("=")[1].replace(';', '').strip()
            components = re.sub(r"[{} ]","",components1).split(",")
            #print("Components : ",components)
            
        searchstatus = re.match(r'\btdb_phases\b', line)
        if searchstatus:
            flg_Phase = 1  
            PHASESline = line
            phases1 = PHASESline.split("=")[1].replace(';', '').strip()
            phases = re.sub(r"[{} ]","",phases1).split(",")
            #print("Phases for calculation : "  , phases)
        
        searchstatus = re.match(r'\bnum_thermo_phases\b', line)
        if searchstatus:
            flg_Num_Phases = 1
            N_PHASESline = line
            Phases1 = N_PHASESline.split("=")[1].replace(';', '').strip()
            N_Phases = int(np.array(re.sub(r"[{} ]","",Phases1).split(","))[0])
            #print("Num of Phases : ",N_Phases)
        
if flg_tdbf == 0:
    print("###########################################")
    print("# Error in Input file: No TDB information #")
    print("###########################################")

if flg_Comps == 0:
    print("##################################################")
    print("# Error in Input file: No components information #")
    print("##################################################")

if flg_Phase == 0:
    print("#############################################")
    print("# Error in Input file: No phase information #")
    print("#############################################")

if flg_Num_Comps == 0:
    print("###########################################")
    print("# Error in Input file: No NUMCOMPONENTS information #")
    print("###########################################")

if flg_Num_Phases == 0:
    print("###########################################")
    print("# Error in Input file: No NUMPHASES information #")
    print("###########################################")

db_sys = Database(tdbfname)

# 'VA' required for pycalphad equilibrium calculations

components = components+['VA']

# In[]:
    
# Checking for presence of components from input file in TDB components 
# component sys contains components from input as present in TDB file
# Usually Components for pycalphad modules are in upper case

component_of_sys = list(db_sys.elements)
component_sys = []

for component in components:
    flg_component = 0
    for comp in component_of_sys:    
        if comp.upper() == component.upper():
            flg_component = 1
            component_sys.append(comp)              
            break
    if flg_component == 0:
        
       print("Error: Component :", component," not component in TDB file")

# In[]:
       
phases_sys = list(db_sys.phases.keys())
phases_sys = sorted(phases_sys) 

print("\n````````````````````````````````")
print("TDB file name : ", tdbfname)
print("Components : ", component_sys)
print("tdb_phases for calculation : ", phases)
print("````````````````````````````````\n")

# In[]:
    
Liq_flg = 0
Sol_flg = 0

for phase in phases_sys:
    
    if phases[0] == phase:    
        Sol_flg = 1
        Solid_Phase = phase
        print('Solid Phase:', Solid_Phase)

    if phases[1] == phase:    
        Liq_flg = 1
        Liquid_Phase = phase
        print('Liquid Phase:', Liquid_Phase)

if Liq_flg == 0:
    
    print('Liquid Phase : ', phases[1],' not a TDB phase')

if Sol_flg == 0:
    
    print('Solid Phase : ', phases[0],' not a TDB phase')
     
# In[]    

# Obtain input Liquid Composition (Alloy Composition) for calculation
# Liquid Phase is last phase in tdb_phases

Liq_id = N_Phases - 1
Liq_id_arr = np.full(N_Phases,Liq_id)

flg_ceq = 0
flg_ceq_Liq = 0

with open('Input.in', 'r') as my_file:    
    
    flines=my_file.readlines()
        
    for line in flines:
 
        searchstatus = re.match(r'\bceq\b', line)
        if searchstatus:
            flg_ceq = 1  
            ceq_line = line
            ceq_tuple = ceq_line.split("=")[1].replace(';', '').strip()
            ceq = np.array(re.sub(r"[{} ]","",ceq_tuple).split(","), dtype = float)
            
            if (ceq[0:N_Phases] == Liq_id_arr).all() :
               
                Liq_Comp = ceq[N_Phases:]
                flg_ceq_Liq = 1
 
if flg_ceq == 0:

    print("###########################################")
    print("# Error in Input file: No ceq information #")
    print("###########################################")
    
elif flg_ceq_Liq == 0:
        
    print("###########################################")
    print("# Error in Input file: No ceq_Liq information #")
    print("###########################################")
    
    
print( 'Liquid composition :', Liq_Comp)    
# In[]:

# comp_sys contains all components except 'VA'    
# Independent components are the components for calculation
# Dependent component is the last components

comp_calc = [x for x in component_sys if x != 'VA']

Comp_Indep = comp_calc[0:-1]  # Component for Calculation
Comp_Depdt = [comp_calc[-1]] # Dependent Component

print("Independent Components : ", Comp_Indep) 
print("Dependent Component : " , Comp_Depdt)

# In[]:
  
def equi_composition_phase(Temperature,Composition,phase,phases_in_sys):  
    
    
    """
     Function to compute equilibrium phase composition at a Temperature 
 
     Inputs:
 
     1. Temperature 
     2. Composition : Alloy Composition in as array (Ex: [0.5], [0.4 0.2])
                      Does not include the dependent component 
     3. phase : phase for composition is found (ex: "LIQUID")  
     4. phases_in_sys : list of phases present in the system to be considered for
                       the calculation (ex: ["FCC_A1","LIQUID"])  
 
     return:
         
     Output array containing equilibrium composition for phase  
     
     Note : Output array doesn't contain composition of dependent component 
            Composition are in the order of elements as taken from .in file
    """
    
    phase_comp = []
    phase = [phase]
    
    conds = {v.T: Temperature, v.P:101325}
    var_comp = [v.X(comp_calc[i]) for i in range(len(comp_calc)-1)]
    conds_comp = {var_comp[i]: Composition[i] for i in range(len(var_comp))}
    
    conds.update(conds_comp)
    
    eq_comp = equilibrium(db_sys, component_sys, phases_in_sys, conds,output='GM')
    Phase_at_equi = np.array(eq_comp.Phase.squeeze())
    Phase_at_equi = [x for x in Phase_at_equi if x] # Eliminate Nan Entries
    Phase_at_equi = sorted(Phase_at_equi)
    
    if phase[0] in Phase_at_equi:
        
        for comps in Comp_Indep:
        
            phase_compo = np.array(eq_comp.X.where(eq_comp.Phase == phase).sel(component = comps).squeeze())
            phase_compo = phase_compo[~(np.isnan(phase_compo))]
            phase_comp.append(phase_compo)
        
        phase_comp = np.array(phase_comp)
        
        return phase_comp 
     
    else:
         
        print('Phase not present for the Temperature and Composition')
        return    

# In[]:
    
#print(equi_composition_phase(800, Liq_Comp, Solid_Phase, phases_sys)[:,0])

# In[]:

def phase_in_equilibrium(Temperature,Composition,phases_in_sys):
    
    """
    
    Function to find phases at equilibrium(sorted) at a Temperature 
    for a given composition 
    
    Input :
    
    1. Temperature
    2. Composition : Alloy Composition in as array (Ex: [0.5], [0.4 0.2])
                   Does not include the dependent component 
    3. phases_in_sys : list of phases present in the system to be considered 
                       for the calculation ((ex: ['HCP_A3','FCC_A1','LIQUID']))  
 
    return :
         
    Output array containing equilibrium phases(sorted)    
    
    """
    
    
    conds = {v.T: Temperature, v.P:101325}
    var_comp = [v.X(comp_calc[i]) for i in range(len(comp_calc)-1)]
    conds_comp = {var_comp[i]: Composition[i] for i in range(len(var_comp))}
    
    conds.update(conds_comp)

    eq_com = equilibrium(db_sys, component_sys,phases_in_sys,conds,output='GM')
    
    Phase_at_equi = np.array(eq_com.Phase.squeeze())
    Phase_at_equi = [x for x in Phase_at_equi if x]
    Phase_at_equi = sorted(Phase_at_equi)
      
    return Phase_at_equi

# In[]:
    
#print(phase_in_equilibrium(734, Liq_Comp, phases_sys)) 

# In[]:

# Find melting point of Liq_Comp    

T_melt, Output_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(400, 2500, Liq_Comp, phases_sys)  
slope = Liq_Tem_ent.sign_slope_Liquidus(Liq_Comp,T_melt,phases)

if slope > 0:

   c_perturb = np.full(N_comps-1,0.0005)

else:

   c_perturb = np.full(N_comps-1,-0.0005)  

del_c = 0.01
  
# In[]:
    
# Liquidus line is specific to a system. Running Liquidus_Temp_entire generates
# a csv file containing Liq Comp and Corresponding Melting Points. 
# if file not present run Liquid_Tem_entire function as done below    

"""  
Liquid_file = df(pnd.read_csv('Liquidus Line Al_Zn',delimiter=","))  
T_Liq = Liquid_file['Temp']
str_name = 'X_'+ Comp_Indep[0] # only for binary
Compo_Liq_arr = Liquid_file[str_name]

LiqT_arr = [ [Compo_Liq_arr[i], T_Liq[i]] for i in range(len(T_Liq)) ]  
LiqT_arr = np.array(LiqT_arr) 
"""
# In[]:

#LiqT_arr has Liq. composition, corresponding Melting Point 
 
LiqT_arr = Liq_Tem_ent.Liquid_Tem_entire(Liq_Comp,del_c,0)
LiqT_arr = np.array(LiqT_arr)

# In[48]:

Compo_Liq_arr = LiqT_arr[:,0].squeeze() # Liq Composition array
T_Liq = LiqT_arr[:,1]                   # Corresponding melting points  
T_Liq = np.array(T_Liq)                   
Compo_Liq_arr =  [float(x[0].item()) for x in Compo_Liq_arr] # for interp1d
T_Liq = LiqT_arr[:,1]                   # Corresponding melting points                   
T_Liq = [float(x) for x in T_Liq]    # for interp1d

# In[]:

# Comp_Sol_arr has Sol. composition in equilibrium with Comp_Liq_arr        
Compo_Sol_arr = np.zeros(len(Compo_Liq_arr))
 
for i in range(len(Compo_Liq_arr)):
 
    Compo_Sol_arr[i] = equi_composition_phase(T_Liq[i],Compo_Liq_arr[i]+c_perturb,Solid_Phase,phases_sys)[:,0]

# In[]:

# Liquidus Temp as a function of Liq Composition
f_cubic_Liquidus = interp1d(Compo_Liq_arr, T_Liq, kind='cubic')

# Solid Composition as a function of Liq Composition 
f_cubic_Liq_Sol = interp1d(Compo_Liq_arr, Compo_Sol_arr, kind='cubic')

# Liquid Composition as a function of Solid Composition 
f_cubic_Sol_Liq = interp1d(Compo_Sol_arr, Compo_Liq_arr, kind='cubic') 

# c_liq_max and c_liq_min are found for the interpolation functions.
# The interpolation is done for liq composition. The corresponding solid 
# composition for the min or max of Liq Composition will be outside the 
# interpolation range for finding melting point      

c_max = np.max(Compo_Liq_arr) - del_c
c_min = np.min(Compo_Liq_arr) + del_c

#print(c_max,c_min)

c_liq_max, c_liq_min = c_max, c_min

c_sol_min = f_cubic_Liq_Sol(c_min) 
c_sol_max = f_cubic_Liq_Sol(c_max) 

if c_sol_max > c_max :
    
    c_liq_max = f_cubic_Sol_Liq(c_max)
    
elif c_sol_min < c_min:
      
    c_liq_min = f_cubic_Sol_Liq(c_min)
       
#print(c_liq_min,c_liq_max)     

# In[]:
    
Ga_mma = GTC.Gibbs_Thom_Coef(350, 1800, Liq_Comp) #m-K Gibbs Thomson Coefficient
print(Ga_mma)
D = 1*10**-9       #m**2/s  # Diffusivity
error_tolerance = 5e-4   # for fsolve errors 

# In[62]:

def solve_eqn(G,C_0,T):
    
    # Find solutions c_l, Pe, x_tip = 1/R_tip
    # Solve using fsolve. 
    
    c_l = G[0]
    
    x_tip = abs(G[1])  # fsolve might throw negative values. 
    Pe = abs(G[2]) 

    if c_l < c_min:        
        c_l = c_liq_min
    
    if c_l > c_max:        
        c_l = c_liq_max
        
    Tm_cl = f_cubic_Liquidus(c_l).item() 
    c_s = f_cubic_Liq_Sol(c_l) 
  #  print(c_s)
    Tm_cs = f_cubic_Liquidus(c_s).item()
    
    d_tip = abs(Ga_mma/(Tm_cs - Tm_cl))*x_tip
    
    eq_1 = abs(2*Ga_mma*x_tip - (Tm_cl - T))
    eq_2 = abs(Pe - 2*(np.pi)**2*d_tip)
    eq_3 = abs((c_l - C_0)/(c_l - c_s) - Pe*(np.exp(Pe))*sp.exp1(Pe))
    
    return [eq_1,eq_2,eq_3]
 
# In[]:    
   
def find_V_R_tip(Guess,C_far,T,num_T_step):
    
    # The solved are better if done in steps for large undercoolings
    
    T_arr = np.linspace(T_melt, T ,num_T_step)
    
    for T_em in T_arr[1:]:
        
        res_T, info, ier, msg = opti.fsolve(solve_eqn,Guess,args = (C_far,T_em),full_output= True)
        res_T = abs(res_T)   # absolute value is the solution
        Guess = res_T
         
    error = solve_eqn(res_T,C_far,T) 
    error = np.array(error)

    R_tip = 1/res_T[1]    
    V_tip = 2*D*res_T[2]*res_T[1]
    
    if error.max() > error_tolerance:

        print('error greater than error_tolerance :', error_tolerance)
        print(error) 
        return V_tip, R_tip, res_T
    
   # print( ' error :' , error)
     
    return V_tip, R_tip, res_T

# In[]:

def VR_entire_range(Guess, C_alloy):

    """    
        Find V tip and R tip for the entire solidification range and write
        to a .csv

        Input : Initial Guess and Alloy Composition
        Melting Point of alloy is found using Temp_Range_equilibrium
       
    """
    
    df_final = df()
    T_melt, Output_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(400, 2500, C_alloy, phases_sys)
    T_Solidus = Output_arr.T_Low.where(Output_arr.Phase_Sol == Solid_Phase)
    print('Initial Guess :', Guess)
    
    T_dif = 0.01 
    T = T_melt - T_dif
    
    del_T = T_melt - T

    V_tip , R_tip, res_1 = find_V_R_tip(Guess,C_alloy,T,2)

    df_Temp = df({'T': T, 'Pe': res_1[2] , 'V_tip': V_tip, 'R_tip' : R_tip, 'c_l' : res_1[0]},index=[0])
    df_final = df_final.append(df_Temp)

    Guess = res_1

    T_start = int(T_melt)
    T_end = int(T_Solidus) + 1

    print('start' , T_start)
    print('End' , T_end)

    T = int(T_start)

    while (T <= T_start and T >= T_end):
    
        T = T - T_dif
        del_T = T_melt - T
        
        if abs(del_T) >= 5:
        
           T_dif = 0.5
            
        V_tip , R_tip, res_T = find_V_R_tip(Guess,C_alloy,T,2)
        
        df_Temp = df({'T': T, 'Pe': res_T[2] , 'V_tip': V_tip, 'R_tip' : R_tip, 'c_l' : res_T[0]},index=[0])
        df_final = df_final.append(df_Temp)
        Guess = res_T
        
    print('Done')    
    
    file_name = 'VR_Tip.txt'
    
    df_final.to_csv(file_name, sep = ' ',index = False)
    
    return df_final

# In[]:

"""    
c_Guess = Liq_Comp[0]
R_Inv_Guess = 10
Pe_Guess = 1e-5
Guess = [c_Guess,R_Inv_Guess ,Pe_Guess]    
VR_entire_range(Guess, Liq_Comp)
"""
# In[]:
"""
del_T = 13
print('Undercooling :' ,del_T)
print('Guess :', Guess)
T = T_melt - del_T
print('Temperature :', T)
V_tip , R_tip, res_1 = find_V_R_tip(Guess,Liq_Comp,T, 400)
print('V_tip :', V_tip)
print('R_tip :', R_tip)   
print(res_1)
"""

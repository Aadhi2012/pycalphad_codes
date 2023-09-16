#!/usr/bin/env python
# coding: 
    
# In[]:

import sys
import numpy as np
from pycalphad import Database, equilibrium
import pycalphad.variables as v
import re
import time 
from operator import itemgetter
from pandas import DataFrame as df

import Liquidus_Tem as Solv_Liq_Temp # To find Melting point 

# In[]:

"""
  
 This script is used to find Liquidus Line for a Solid Liquid Equilibrium
 when melting point for one composition is known. written only for binary 
 system.
 
 Using secant method the melting point for a new composition is found. (found
 for composition close to composition for which melting point is known. 
 Will work for any composition if bounds for temperature is known. 
 The bounds would be Melting point of component which is higher than other 
 component and Eutectic temp. Here no bounds are included)
  
 Input from User for the script:
     
          1: Database file(.tdb) 
          2: Solid Phase in equilibrium with liquid at the melting point 
          3: Liquid Phase
          4: Components for Calculation 
          5: Alloy Composition
     
"""
     
# In[]:
    
"""    

Inputs read from .in input file taken as argument using sys.argv

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
            print("Tdb file name : ", tdbfname)

        searchstatus = re.match(r'\bNUMCOMPONENTS\b', line)
        if searchstatus:
            flg_Num_Comps = 1
            N_COMPONENTSline = line
            compts = N_COMPONENTSline.split("=")[1].replace(';', '').strip()
            N_comps = int(np.array(re.sub(r"[{} ]","",compts).split(","))[0])
            print("Num of Comps : ",N_comps)
        
        searchstatus = re.match(r'\bCOMPONENTS\b', line)
        if searchstatus:
            flg_Comps = 1
            COMPONENTSline = line
            components1 = COMPONENTSline.split("=")[1].replace(';', '').strip()
            components = re.sub(r"[{} ]","",components1).split(",")
            print("Components : ",components)
            
        searchstatus = re.match(r'\btdb_phases\b', line)
        if searchstatus:
            flg_Phase = 1  
            PHASESline = line
            phases1 = PHASESline.split("=")[1].replace(';', '').strip()
            phases = re.sub(r"[{} ]","",phases1).split(",")
            print("Phases for calculation : "  , phases)
        
        searchstatus = re.match(r'\bnum_thermo_phases\b', line)
        if searchstatus:
            flg_Num_Phases = 1
            N_PHASESline = line
            Phases1 = N_PHASESline.split("=")[1].replace(';', '').strip()
            N_Phases = int(np.array(re.sub(r"[{} ]","",Phases1).split(","))[0])
            print("Num of Phases : ",N_Phases)
        
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
  #      print('Solid Phase:', Solid_Phase)

    if phases[1] == phase:    
        Liq_flg = 1
        Liquid_Phase = phase
   #     print('Liquid Phase:', Liquid_Phase)

if Liq_flg == 0:
    
    print('Liquid Phase : ', phases[1],' not a TDB phase')

if Sol_flg == 0:
    
    print('Solid Phase : ', phases[0],' not a TDB phase')
     
# In[]    

# Obtain input Liquid Composition (Alloy Composition) for calculation
# Liquid Phase is last phase in tdb_phases so Liq_id is N_Phases - 1

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

Temp_Tol = 1e-4      # Temperature tolerance for searches
Comp_Tol = 5e-5     # composition tolerance
T_step = 1e-4        # For finding slope 
    
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
       the calculation (ex: ["FCC_A1","LIQUID"], ["FCC_A1","LIQUID", "HCP_A3"])
 
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
    
#print(equi_composition_phase(699.2398, [0.21], Liquid_Phase, phases_sys)[:,0])

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
       for the calculation(ex: ["FCC_A1","LIQUID"], ["FCC_A1","LIQUID", "HCP_A3"]) 
 
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
    
#print(phase_in_equilibrium(699, [0.21], phases_sys)) 

# In[]:

def sign_slope_Liquidus(C,T,phases_calc_sort):
    
    """
    
    Function to find Liquidus slope for finding c_perturb
    
    Inputs:
    
        C : Composition to find Liquidus line slope at Temperature T 
        At ( C,T ) Solid + Liquid equilibria should exist
        phases_calc_sort : phases in Equilibrium 
    
        c_perturb: For pycalphad equilibrium calculations 
                   A point on Liquidus line has equilibrium phase LIQUID not
                   LIQUID + SOLID,hence to find Composition of SOLID, 
                   c_perturb is found so that equilibrium phasbe will be 
                   LIQUID + SOLID for C + c_perturb
    
    """
    
    global c_perturb
    
    T_for = T - T_step
    c_calc = C
    phases_calc_sort = sorted(phases_calc_sort)
    phases_C = phase_in_equilibrium(T,c_calc,phases_sys)
    #print(phases_C)
    
    if phases_C == phases_calc_sort:
    
        Liq_compo_1 = equi_composition_phase(T,c_calc,Liquid_Phase,phases_sys)[:,0]
        Liq_compo_2 = equi_composition_phase(T_for,c_calc,Liquid_Phase,phases_sys)[:,0]    
    
        Liq_slope = (Liq_compo_1 - Liq_compo_2)/T_step
        
        if Liq_slope > 0:

            c_perturb = np.full(N_comps-1,0.0005)

        else:

            c_perturb = np.full(N_comps-1,-0.0005)  
            
        return Liq_slope
    
    else:
        
        print('Equilibrium Not Found')
        
        return 0

# In[]:

def x_Secant(c_alloy,c_slope,T): 
    
     
    """
    Function to find new value using Secant method
    
    Inputs:
    
    c_slope : composition to compute slope. Should have Liq + Solid equilibria
    T : Liquidus Temperature of c_slope
    c_perturb : pertubation to find composition of phases in equilibrium 
    c_alloy : alloy composition to solve (c_liq - c_alloy) = 0
    
    Outputs:
    
    Liquidus Temperature for c_alloy    
    
    """
    
    T_for = T - T_step  
    
    Liq_compo_1=equi_composition_phase(T,c_slope+c_perturb,Liquid_Phase,phases_sys)[:,0]
    Liq_compo_2=equi_composition_phase(T_for,c_slope+c_perturb,Liquid_Phase,phases_sys)[:,0]    
    
    slope_Liq = T_step/(Liq_compo_1 - Liq_compo_2)
    
    x = float(new_x(slope_Liq, T, Liq_compo_1 - c_alloy))
    
    return x

# In[]:

def new_x(slope,x_known,f_known):
    
    """
    x is found using Secant method
    
    """
    x = x_known - f_known*(slope)
    
    return x 

# In[]:

def find_c_equilibrium(c_alloy,C_eval,c_sol_eval,T_new,phases_calc_sort):
    
    """
    Function to find composition for equilibrium calculation at a Temperature
    This is required for pycalphad equilibrium calculations in Solid +
    Liquid region
    
    c_alloy : Alloy Composition
    c_sol_eval : half of Composition of solid + liquid.
    C_eval: composition for which Liquidus temperature to be found
    T_new: Temperature at which the composition for equilibrium calculation 
           has to be found
    phases_calc_sort:  [Solid,Liquid] in equilibrum sorted
    
    """
    
    ph_alloy = phase_in_equilibrium(T_new,c_alloy, phases_sys) 
    ph_eval = phase_in_equilibrium(T_new,C_eval, phases_sys)
    ph_sol = phase_in_equilibrium(T_new,c_sol_eval, phases_sys)
    
    Solid_Phase =  [x for x in phases_calc_sort if x!= Liquid_Phase][0]
    
    ph_all = [ph_alloy,ph_eval,ph_sol]
    comp_all = [c_alloy,C_eval,c_sol_eval]
    
    C_S = np.full(N_comps-1,-1)
    C_L = np.full(N_comps-1,-1)
    
    for i in range(len(ph_all)):
        
        if ph_all[i] == phases_calc_sort:
              
            c = comp_all[i]
        
            return c
    
    #  If no equilibria found for above three compositions. C_S and C_L are 
    #  found to do a binary search between C_S and C_L
    
    for i in range(len(ph_all)):
        
        if ph_all[i] == [Solid_Phase]:
    
            C_S = comp_all[i]
            
        
        if ph_all[i] == [Liquid_Phase]:
    
            C_L = comp_all[i]
            
    # Binary Search when C_S and C_L are found
     
    if (C_S != np.full(N_comps-1,-1) and C_L != np.full(N_comps-1,-1)):
        
        c = Compo_Binary_Search(C_S,C_L,T_new,phases_calc_sort)
        return c
        
    else:
        
        print('Composition for equilibrium at ',T_new,' not found')
        
    return -1 

# In[]:

def Compo_Binary_Search(C_S,C_L,Temp,phases_calc_sort):
    
    """ 
    Function to do a Binary Search between C_S and C_L at a Temperature to 
    find Solid + Liq equilibrium 
   
    Input: 
   
    C_S: phase at equilibrium : Solid
    C_L: phase at equilibrium : Liquid
     
    Temp : Temp at which Composition has to be found
    phases_calc_sort:  [Solid,Liquid] in equilibrum sorted
   
    return:
        
    composition      
    
    """
    
    c_cur = (C_S + C_L)*0.5
    phase_cur = phase_in_equilibrium(Temp, c_cur, phases_sys) 
    
    while phase_cur != phases_calc_sort:
                
        if phase_cur == [Liquid_Phase] :
    
            C_L = c_cur
    
        else:
             
            C_S = c_cur
        
        if (np.max(abs(C_S - C_L)) < Comp_Tol):
            
            print('Composition for Equilibrium at ',Temp ,'not found')
            return -1
        
        c_cur = (C_S + C_L)*0.5
        phase_cur = phase_in_equilibrium(Temp, c_cur, phases_sys)       
    
    return c_cur

# In[]:

def find_Liquidus_Temperature(c_alloy,c_known,T_known,phases_calc_sort):
    
    """
    
    Function to find Liquidus temperature for c_alloy 
    
    Inputs:
    
    c_alloy: Alloy Composition for which Liquidus Temperature has to be found 
    c_known: composition for which Liquidus Temperature(T_known) is known 
    
    phases_calc_sort:  [Solid,Liquid] in equilibrum sorted
    
    return:
        Liquidus temperature and phases in equilibrium
        
    Liquidus Temperature , [Solid, Liquid] sorted     
    
    """
 
    C_eval = c_known
    T_eval = T_known
    
    # c_sol_eval is the (solid composition + Liquid composition)* 0.5
    # At c_sol_eval Solid + Liquid exists for nearby temperatures above and 
    # below

    c_sol_eval = 0.5*(equi_composition_phase(T_eval,c_known+c_perturb, Solid_Phase, phases_sys)[:,0]+ C_eval)
    comp_diff = c_alloy - c_known
    
    if (np.max(abs(comp_diff))) <= Comp_Tol:
                
        T_new = T_eval 
        
    # Melting point of c_alloy is found using secant method. New value of Temp
    # after every iteration will have a Liquid Composition which will be 
    # matched with c_alloy within Comp_Tol     
        
    while np.max(abs(comp_diff)) > Comp_Tol :   
        
        T_new = float(x_Secant(c_alloy,C_eval,T_eval))
        C_eql = find_c_equilibrium(c_alloy,c_known,c_sol_eval,T_new,phases_calc_sort)
        # At C_eql Liquid + SOlid for T_new exists 
        
        if C_eql != -1:
        
            c_liq_new = equi_composition_phase(T_new,C_eql,Liquid_Phase, phases_sys)[:,0]
            # C_liq_new Liq Comp at T_new
            comp_diff = abs(c_liq_new - c_alloy)
            T_eval = T_new
            C_eval = c_liq_new
        
        else:
            
            break
        
    phase_equi_cur = phase_in_equilibrium(T_new,c_alloy+ c_perturb, phases_sys)
    
    if phase_equi_cur == phases_calc_sort:
            
        #print('Liquidus Temperature Found : ' ,T_new)     
        return T_new,phase_equi_cur
        
    else:
           
        print('Liquidus Temperature not Found for : ', c_alloy)             
        
        return T_new,phase_equi_cur 

# In[]:

#T_melt, Output_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(400, 2500, Liq_Comp, phases_sys)
#slope = sign_slope_Liquidus(Liq_Comp,T_melt,phases)

#print(T_melt, Output_arr)
#print(c_perturb)

# In[]:

#T_fin = find_Liquidus_Temperature([0.259437],Liq_Comp,T_melt,phases) 
#print(T_fin)

# In[]:
    
def Liquid_Tem_entire(C_alloy,del_c, flag_write_file):
    
    """
      The function finds writes composition and corresponding liquidus 
      temperature for Liquidus line using find_Liquidus_Temperature() 
      using known liquidus temperature
      
      Inputs:
      
          C_alloy: Alloy Composition for which Liquidus Temp is known 
          del_c = composition step at which Liquidus Temp has to be found
          
          flag_write_file : If true writes Liquidus Temp and Comp onto a csv
          Since Liquidus line is specific for a system from the phase diagram
          it can be written onto a file to use for other calculations.
          
      return:
          
          Array containing [Composition, Liquidus Temp]
          Composition is type array  
          Output : Array[Composition, Liquidus Temp] over entire range 
          
    """
    
    del_c_arr = np.full( N_comps-1 ,del_c)
    tim = time.time()    
                
    T_melt, Output_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(400, 2500, C_alloy, phases_sys)
    slope = sign_slope_Liquidus(C_alloy,T_melt,phases)
    # slope is called to initialize c_perturb
    
    Compo_LiqT_arr = [] 
    Compo_LiqT_arr.append([C_alloy,T_melt])
    
    c_know = C_alloy
    T_know = T_melt
    phase_in_equil = phases
    
    # The composition are evaluated after rounding of compositions to x th
    # digit after decimal.  x determined by del_c
    
    round_digit = int(np.log10(int(1/del_c)))
    round_c_alloy = [ round(c,round_digit) for c in C_alloy]
    
    c_round = C_alloy
    T_Liq_round = T_melt
    
    if abs(C_alloy - round_c_alloy) > 1e-5:
        
        c_round = round_c_alloy
        T_Liq_round, phase_in_equil = find_Liquidus_Temperature(c_round,c_know,T_know,phases)
        c_know = c_round    
        T_know = T_Liq_round 
        Compo_LiqT_arr.append([c_round,T_Liq_round])
    
    print('Increasing Composition in steps of ',del_c)
    
    while (phase_in_equil == phases and np.max(c_know) < 0.98):
     
        c_find = c_know + del_c_arr     
        T_Liq_res, phase_in_equil = find_Liquidus_Temperature(c_find,c_know,T_know,phases)
         
        if phase_in_equil == phases:
        
           Compo_LiqT_arr.append([c_find,T_Liq_res])
        
        c_know = c_find
        T_know = T_Liq_res
        
    c_know = round_c_alloy    
    T_know = T_Liq_round   
    
    print('Decreasing Composition in steps of ',del_c)  
    
    del_c_arr = -del_c_arr
    phase_in_equil = phases
    
    while (phase_in_equil == phases and np.min(c_know) > 0.02) :
     
        c_find = c_know + del_c_arr   
        T_Liq_res, phase_in_equil=find_Liquidus_Temperature(c_find,c_know,T_know,phases)
           
        if phase_in_equil == phases:
        
           Compo_LiqT_arr.append([c_find,T_Liq_res])
        
        c_know = c_find
        T_know = T_Liq_res
   
    print('Time taken : ',time.time()- tim,' --- s') 
    
    res_arr = Compo_LiqT_arr
    LiqT_arr = sorted(res_arr, key=itemgetter(0)) # sort according to 1st element
    LiqT_arr = np.array(LiqT_arr)
    
    if flag_write_file == 1:
        
        Compo_Liq_arr = LiqT_arr[:,0].squeeze() # Liq Composition array
        T_Liq = LiqT_arr[:,1]                   # Corresponding melting points                   
        Compo_Liq_arr = [x[0].item() for x in Compo_Liq_arr] # only for binary
        str_name = 'X_'+ Comp_Indep[0] # only for binary
        df_Melts = df({str_name: Compo_Liq_arr, 'Temp': T_Liq})       
        df_Melts.to_csv('Liquidus Line Al_Zn', index= False)

    return LiqT_arr
            
# In[]:    

#res_arr = Liquid_Tem_entire(Liq_Comp,0.01,1)    
    

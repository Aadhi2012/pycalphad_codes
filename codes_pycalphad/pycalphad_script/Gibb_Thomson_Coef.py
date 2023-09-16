#!/usr/bin/env python
# coding: utf-8

# In[]:

import sys
import numpy as np
import sympy as sym
from pycalphad import Model, Database, equilibrium
import pycalphad.variables as v
import re

import Liquidus_Tem as Solv_Liq_Temp

# In[]:

"""

 Script to evaluate Gibbs-Thomson Coefficient. 
 
 Input for Temp_Range_Equilibrium: 
    
    Alloy Composition, phases and Solid phase in equilibrium with Liquid Phase
  
 Input from User for the script:
     
          1: Database file(.tdb) 
          2: Solid Phase in equilibrium with liquid at the melting point 
          3: Liquid Phase
          4: Components for Calculation 
          5: Alloy Composition
          6: Surface Tension          # In SI units (J/m2)
 
 This .py Script uses Temp_Range_Equilibrium() from Liquidus_Tem
 imported as Solv_Liq_Temp for finding melting point for a binary alloy  
         
"""

# In[]: 

"""    

Inputs read from .in input file taken as argument using sys.argv

Inputs read from input file:
    
tdbfname, NUMCOMPONENTS, COMPONENTS, tdb_phases, num_thermo_phases, ceq, V, GAMMA    

"""
print('Reading Input file')

#with open(sys.argv[1], 'r') as my_file:    
with open('Input.in', 'r') as my_file:    
    
    flines=my_file.readlines()

    flg_tdbf, flg_Comps, flg_Num_Comps, flg_Phase, flg_Num_Phases = 0,0,0,0,0  
    flg_surf_tens, flg_V = 0, 0
    
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
        
        searchstatus = re.match(r'\bGAMMA\b', line)
        if searchstatus:
            flg_surf_tens = 1
            Gammaline = line
            Gamma1 = Gammaline.split("=")[1].replace(';', '').strip()
            GAMMA = np.array(re.sub(r"[{} ]","",Gamma1).split(","),dtype = float)[0]
            print('Surface Tension :', GAMMA)
            
        searchstatus = re.match(r'\bV\b', line)
        if searchstatus:
            flg_V = 1
            Vline = line
            V_m = Vline.split("=")[1].replace(';', '').strip()
            V_m = float(V_m)
            print('Molar volume : ', V_m)

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

if flg_surf_tens == 0:
    print("###########################################")
    print("# Error in Input file: No GAMMA information #")
    print("###########################################")

if flg_V == 0:
    print("###########################################")
    print("# Error in Input file: No molar value information #")
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
 #       print('Solid Phase:', Solid_Phase)

    if phases[1] == phase:    
        Liq_flg = 1
        Liquid_Phase = phase
  #      print('Liquid Phase:', Liquid_Phase)

if Liq_flg == 0:
    
    print('Liquid Phase : ', phases[1],' not a TDB phase')

if Sol_flg == 0:
    
    print('Solid Phase : ', phases[0],' not a TDB phase')
     
# In[]    

# Obtain input Liquid Composition (Alloy Composition) for calculation

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

# comp_calc contains all components except 'VA'    
# Independent components are the components for calculation
# Dependent component is the last components

comp_calc = [x for x in component_sys if x != 'VA']

Comp_Indep = comp_calc[0:-1]  # Component for Calculation
Comp_Depdt = [comp_calc[-1]] # Dependent Component

print( 'Independent Components :', Comp_Indep)
print( 'Dependent Component :', Comp_Depdt)    

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
    
#print(equi_composition_phase(734, Liq_Comp, Liquid_Phase, phases_sys))

# In[]:

def phase_in_equilibrium(Temperature,Composition,phases_in_sys):
    
    """
    
    Function to find phases at equilibrium(sorted) at a Temperature 
    for a given composition 
    
    Input :
    
    1. Temperature
    2. Composition : Alloy Composition in as array (Ex: [0.5], [0.4 0.2])
                   Does not include the dependent component 
    3. phases_in_sys : list of phases present in the system to be considered for
       the calculation (ex: ["FCC_A1","LIQUID"], ["FCC_A1","LIQUID", "HCP_A3"]) 
 
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

def der_GM_T(Phase,comps,Temp):
    
    """
    
    Function returns value of derivative of free energy(G) of phase
    w.r.t T evaluated at comps 
    
    Input:
        
        Phase : phase for which derivative is evaluated
        comps : composition w.r.t elements present in order as in input file
        Temp : Temp
        
    return:    
    
        returns the value of the derivative    
    
    """
    
    phase_mod = Model(db_sys, component_sys, Phase)
    my_list = [v.Y(Phase,0,comp_calc[i]) for i in range(len(comp_calc))]
    my_list.append(v.T)

    GM_phase = phase_mod.GM
    GM_der = GM_phase.diff(my_list[-1])    
    
    subs_arr = [(my_list[i],comps[i]) for i in range(len(my_list)-2)] 
    subs_arr.append((my_list[-2] , 1 - np.sum(comps)))         
    subs_arr.append((my_list[-1] , Temp))
    
    der_val = GM_der.subs(subs_arr)
    der_val = sym.simplify(der_val)
    
    return der_val 

# In[]:

def Gibbs_Thom_Coef(start, stop, C_alloy):
    
    """
    Function to find Gibbs Thomson Coefficient units : mol-K/m**2
    
    Inputs:
           C_alloy : Alloy Composition
           
    return:
        
    Value of Gibbs Thomson Coefficient in K-m 
    
    """
    
    T_melt, out_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(start,stop, C_alloy,phases_sys)
    
    Solids_at_equi = np.array(out_arr.Phase_Sol)
    
    if Solid_Phase in Solids_at_equi:
    
       Liq_Com = equi_composition_phase(T_melt,C_alloy,Liquid_Phase,phases_sys)[:,0]
       Sol_Com = equi_composition_phase(T_melt,C_alloy,Solid_Phase,phases_sys)[:,0]
       
       val_liq = der_GM_T(Liquid_Phase, Liq_Com, T_melt)
       val_sol = der_GM_T(Solid_Phase, Sol_Com,  T_melt)
    
       val_dif = val_sol - val_liq
    
       Gib_Thom_coef = GAMMA/ val_dif*V_m
   
    return Gib_Thom_coef

# In[]:

#print(Gibbs_Thom_Coef(350,2500,Liq_Comp))
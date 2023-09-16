#!/usr/bin/env python
# coding: 
    
# In[]:

"""

 Script to find the melting point of the alloy and solidication range
 (Solid Phase + Liquid Phase) of a given alloy composition 
 for a two-component system. 
         
 INPUT PARAMETERS for the script    
  
 1. Database File
 2. Solid Phase ( Example : FCC_A1 or BCC_A2  or HCP_A3)
 3. Liquid Phase ( Example :  LIQUID )  
 4. Components for Calculation (component w.r.t which alloy composition is given)
 5. Alloy Composition.
 
"""

# In[]:

import sys
import numpy as np
from pycalphad import Database, equilibrium
import pycalphad.variables as v
from pandas import DataFrame as df
import re

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

print( 'Independent Components :', Comp_Indep)
print( 'Dependent Component :', Comp_Depdt)

# In[]:

Temp_Tol = 1e-4      # Temperature tolerance for Binary Search
    
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
    
#print(equi_composition_phase(830, Liq_Comp, Solid_Phase, phases_sys))

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
                       for the calculation. 
                      (ex: ["FCC_A1","LIQUID"], ["FCC_A1","LIQUID", "HCP_A3"]) 
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
    
#print(phase_in_equilibrium(830, Liq_Comp, phases_sys)) 

# In[]:

def Last_Occurence_Phase(T_Low, T_High, c_alloy, phase_check, phases_in_range): 
       
    """
    
    Binary Search to find Temperature when last phase(taken as phase_check)
    exist when cooled 
    
    Inputs:
        
       1. T_Low : Temperature at which phase_check doesn't exist                                                     
       2. T_High : Temperature at which only phase_check phase exist 
       3. c_alloy : Alloy Composition in as array (Ex: [0.5], [0.4 0.2])
                   Does not include the dependent component 
       4. phase_check : Phase for which last occurence is checked 
                        Ex: ['LIQUID'], ['FCC_A1','LIQUID']                        
       5. phases_in_range : list of phases present in the system to be 
          considered for the calculation.  
         (ex: ["FCC_A1","LIQUID"], ["FCC_A1","LIQUID", "HCP_A3"]) 
    
    Search stops when Temp_diff acheives tolerance and Temp for equilibrium
    is found(flag = 1)  
    
    return:
        
        Function returns an array containing following entries:
        
        Array : [T_Low, T_High, flag]
        
    """
    
    flag = 0
    Temp_diff = abs(T_High - T_Low)
    
    T_Current = (T_High + T_Low)*0.5
    
    while(not(bool(abs(Temp_diff) < Temp_Tol and flag == 1))):  
        
        phase_equi_cur=phase_in_equilibrium(T_Current,c_alloy,phases_in_range)
        
        if phase_equi_cur == phase_check:
            
            flag = 1
            T_High = T_Current
        
        else:
        
            flag = 0
            T_Low = T_Current
        
        T_Current = (T_High + T_Low)*0.5
        Temp_diff = abs(T_High - T_Low)     
    
    if(flag == 0 and Temp_diff < Temp_Tol):
        
        print('Temperature for Last occurence for phase not found' )
        print(' ' )
        print('Try diff start , stop Temperatures' )
        
    return [T_Low,T_High,flag]                

# In[]:

def Temp_Range_Equilibrium(T_Low, T_High, c_alloy, phases_in_range):
    
    """
    
    Function to find Temp range over which Liq + Solid equilibria exists
    Liquidus and Solidus Temp for an alloy composition is found for each 
    (Solid + Liquid) equilibria that exists in the solidification range
    
    Inputs:
    
    1. T_Low : Temperature at which Liquid doesn't exist (Low Temp)
    2. T_High  : Temperature above which only Liquid exist (High Temp)
    3. c_alloy : Alloy Composition in as array (Ex: [0.5], [0.4 0.2])
                   Does not include the dependent component 
    4. phases_in_range : list of phases present in the system 
       to be considered for the calculation.
       (ex: ["FCC_A1","LIQUID"], ["FCC_A1","LIQUID", "HCP_A3"]) 
    
    return:
        
        Returns range of temperature for each equilibria
        Output : Melt_Temp, out_put (x_array)
        
    """
    
    out_put = df()
    df_temp = df()
    i = 0
    flg_equilibrium = 0
    
    # Finding melting point (Liquidus Temperature)
    Liq_Temp_res = Last_Occurence_Phase(T_Low, T_High, c_alloy, [Liquid_Phase],phases_in_range) 
    
    if Liq_Temp_res[-1] == 0:
    
        return
    
    Melt_Temp =  Liq_Temp_res[0]
    #print('Melting Point :', Melt_Temp)
    
    phase_calc_sort = phase_in_equilibrium(Liq_Temp_res[0], c_alloy, phases_in_range)
    res_1 = Liq_Temp_res
    
    #Finding range of Temp for which Equilibria phases exist  
    while sorted(phases) == phase_calc_sort:
    
    # To find all solids in equilibria with Liquid and corresponding temp range    
    # of equilibrium change check statement in while as while Liquid in phases_calc_sort 
    
          flg_equilibrium = 1
        
          Sol_cur = Solid_Phase
          res = Last_Occurence_Phase(T_Low, res_1[0], c_alloy, phase_calc_sort, phases_in_range)
          
          df_temp = df({'Phase_Sol': Sol_cur, 'Phase_Liq': Liquid_Phase, 'T_Low': res[1], 'T_High' :res_1[0]}, index=[i])
          out_put = out_put.append(df_temp)
          i +=1       
          phase_calc_sort = phase_in_equilibrium(res[0], c_alloy, phases_in_range)
          res_1 = res   
    
    if flg_equilibrium == 0:
        
        print("###########################################")
        print("# Error : Equilibrium : ", phases    ," not found #")
        print("###########################################")

    output_arr = out_put.to_xarray()
            
    return Melt_Temp, output_arr

# In[]:

# change T_Low and T_High according melting point of components
# Generally can be set to T_Low = 350 and T_High = 2500 for most systems     

#T_melt, Output_arr = Temp_Range_Equilibrium(350, 2500, Liq_Comp, phases_sys)
#print(T_melt, Output_arr)

# In[]:
    
#Solids_at_equi = np.array(Output_arr.Phase_Sol.squeeze())
#Solids_at_equi = [x for x in Solids_at_equi if x] # Eliminate Nan Entries     
#print(Solids_at_equi)
#print(Output_arr.T_High.where(Output_arr.Phase_Sol == Solid_Phase))    

#!/usr/bin/env python
# coding: utf-8


# In[]:

"""

 Script to write equilibrium compositions and Hessians of free Energy 
 to a .txt file for each phase in tdb_phases from .in file
  
 tdb_phases : ex : ['FCC_A1', 'LIQUID']

 tdb_phases contains phases for whose highest temperature of existence 
 equilibrium is required for writing hessians and equilibrium compositions.  
 
 The Melting point of alloy is required(Temperature at which first Liquid forms)
 
 For binary alloy systems, this .py Script uses Temp_Range_Equilibrium() 
 from Liquidus_Tem module imported as Solv_Liq_Temp   
 
 INPUT PARAMETERS for the script    
  
 1. Database File
 2. Solid Phase ( Example : FCC_A1 or BCC_A2  or HCP_A3)
 3. Liquid Phase ( Example :  LIQUID )  
 4. Components for Calculation (components w.r.t which alloy composition is given)
 5. Alloy Composition.
      
"""
# In[]:
 
import sys    
import numpy as np
from pandas import DataFrame 
from pycalphad import Model, Database, equilibrium
import pycalphad.variables as v
import sympy as sym
from sympy.utilities.codegen import codegen
import re

import Liquidus_Tem as Solv_Liq_Temp 
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
    
#X = equi_composition_phase(850, Liq_Comp, Liquid_Phase, phases_sys)[:,0]
#print(X)

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
    
#print(phase_in_equilibrium(600, Liq_Comp, phases_sys)) 

# In[]:

def Hess_eq(phase):
    
    """ 
    
    Function that returns Hessian of free energy function in symbolic form
    
    Input : 
        
        phase : phase for which hessian is evaluated 
    
    return : return variable : Hess_phase 
        
        Symbolic form of hessian.
        Hessian has dimension (N_comps-1)X(N_comps-1)
        
        Hess_phase[i][j] is function d (dG/d(x_i)) /d(x_j)    
    
    """
    
    phase_mod = Model(db_sys, component_sys, phase)
    
    var = [v.Y(phase, 0, comp_calc[i]) for i in range(len(comp_calc))]
    var.append(v.T)
    
    # var[-1] is 'T'
    # var[-2] is Dependent Component which is var[len(comp_calc)-1]
    
    GM_phase = phase_mod.GM
    Mu = {}
    Hess_phase = {}
    
    for n in var[0:len(comp_calc)-1]:
       
        Mu[n] = GM_phase.diff(n) - GM_phase.diff(var[-2])
    
    for i in var[0:len(comp_calc)-1]:
        
        Hess_phase[i] = {}
        for j in var[0:len(comp_calc)-1]:
            
            Hess_phase[i][j] = Mu[i].diff(j) - Mu[i].diff(var[-2]) 

    return Hess_phase 

# In[]:

#dg_che = Hess_eq(Liquid_Phase)

# In[]:

def compute_hess(Temp, Composition, phase):
    
    """
    Compute hessian for the phase having Composition at Temp
    Hess_all_phases: list contains Hessians for all phases in tdb_phases     
    
    Composition : Alloy composition w.r.t independent components
    
    Output : Dictionary containing hessian for the phase evaluated at Temp for 
             Composition in the order below
    order: 11 , 22, 33, 12, 13, 23 if Independents components were 1,2,3    
    
    """
    
    hess_phase = Hess_all_phases[phase]
          
    val_dict = dict()
    val_dict.update({'T' : Temp})
    
    my_list = [v.Y(phase,0,comp_calc[i]) for i in range(len(comp_calc))]
    my_list.append(v.T)
    
    subs_arr = [(my_list[i],Composition[i]) for i in range(len(my_list)-2)] 
    subs_arr.append((my_list[-2] , 1 - np.sum(Composition))) # Dependent comp      
    subs_arr.append((my_list[-1] , Temp))
    
    for ind_i, i in enumerate(my_list[0:len(comp_calc)-1]):
        
        Hess_val = hess_phase[i][i].subs(subs_arr)
        str_name = 'HSN('+comp_calc[ind_i]+','+comp_calc[ind_i]+')@'+phase
        val_dict.update({str_name:[Hess_val]})   
     
    for ind_i,i in enumerate(my_list[0:len(comp_calc)-2]):
        for ind_j,j in enumerate(my_list[i+1:len(comp_calc)-1]):
            
            str_name = 'HSN('+comp_calc[ind_i]+','+comp_calc[ind_i]+')@'+phase
            Hess_val = hess_phase[j][i].subs(subs_arr)
            val_dict.update({str_name:[Hess_val]})   
    
    #print(val_dict)
    
    return val_dict

# In[]:

def Tem_equi_composition(T_High, c_alloy):
    
    """
    Function to write equi.compo(in the order as in .in file) and 
    Hessians of each phase in phases evaluated at equi.compositions 
    at the Temperature 
    
    Equilibrium Composition written in a single file
    for each phase in phases a Hessian file is written
    
    Inputs :
        
      T_High : temperature at which equilibria starts for the alloy composition 
      c_alloy: alloy composition
      
    return:   
      
      Equilibrium Composition dataframes Containing outputs in order as in 
      phases
        
    """
 
    df_phase_comp = DataFrame() # df that will be written as .txt file
    
    Temp = T_High
    Phase_equili_T = phase_in_equilibrium(Temp, c_alloy, phases_sys)
    
    print('Starting')
    print('Start Temp : ', Temp)
    print('Phases in equilibrium : ', Phase_equili_T)
    # Hessian and Equilibrium composition written till phases at equilibrium is
    # same as phases(sorted)
    
    while Phase_equili_T == phases_sort:
    
          Comps_T_df = DataFrame() # df for comps at that temp
          
          dict_Temp = dict()  # dictionary that'll be converted to Comps_T_df
          dict_Temp.update({'T' : Temp})
          Hess_Temp = {} # List of dictionaries for Hessian from compute_hess
          
          for phase in phases:    
          
              phase_comp = equi_composition_phase(Temp, c_alloy, phase, phases_sys)
              phase_name = phase
              comp_phase = phase_comp[:,0]
              comp_name = ['X_'+ comp_calc[i]+'_' + \
                               phase_name for i in range(len(comp_calc)-1)]
              comp_dict = {comp_name[i]: [comp_phase[i]] for i in range(len(comp_name))}                
              dict_Temp.update(comp_dict)
          # comp_dict has component : composition in the order present in .in file  
              Hess_Temp[phase_name] = compute_hess(Temp, comp_phase, phase)                  
              
          Comps_T_df  = DataFrame.from_dict(dict_Temp)
          df_phase_comp = df_phase_comp.append(Comps_T_df)       
          phases_Temp =  list(Hess_Temp.keys())
        
          for ph_ase in phases_Temp: 
          
              df_cur = (DataFrame.from_dict(Hess_Temp[ph_ase])) 
              # df_dict declared globally
              df_dict[ph_ase] = df_dict[ph_ase].append(df_cur)
                           
          Temp = Temp - 1     
          Phase_equili_T = phase_in_equilibrium(Temp, c_alloy, phases_sys)
    
    for ph_ase in phases_Temp: 
    
        file_str = 'Hess_' + ph_ase + '.txt'  
        df_dict[ph_ase].to_csv(file_str, sep = ',',index = False)
        
    str_file ="Eqb_Comps.txt"
    df_phase_comp.to_csv(str_file, sep = ',',index = False)
    
    print('done')    
    return df_phase_comp

# In[]:

#df_names contain names of hessian dataframes
#df_list contain empty hessian dataframes
df_names = [  x   for x in phases]
df_list = [DataFrame() for x in phases]

# df_dict : dictionary of (df_name : df_list)
df_dict = dict(zip(df_names, df_list))
#print(df_dict)

# In[]:
    
Melt_Temp, res_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(350,1800,Liq_Comp,phases_sys)

print('Melting Point : ' , Melt_Temp)

# In[]:

phases_sort = sorted(phases)    # Phases Sorted       

# T_High : Highest temperature for which equilibrium (phases_sort) exists 
    
T_High = np.array(res_arr.T_High.where(res_arr.Phase_Sol == Solid_Phase))
T_High = T_High[~(np.isnan(T_High))]
T_High = int(T_High)

T_Low = np.array(res_arr.T_Low.where(res_arr.Phase_Sol == Solid_Phase))
T_Low = T_Low[~(np.isnan(T_Low))] 

if (T_Low - int(T_Low)) != 0:

    T_Low = int(T_Low) + 1

T_Low = int(T_Low) 

# In[]:
    
Hess_all_phases = {}

# Hess_all_phases : list of symbolic forms of hessian

for phase in phases:
    
    Hess_all_phases[phase] = Hess_eq(phase)

# In[]:

df_comps = Tem_equi_composition(T_High, Liq_Comp)

# In[]:

def gen_c_phase_eq(phase_equil):
     
    """
    
    Function to generate c codes with Free Energy(G), Mu and Hessian for 
    each phase in phase_equil     
    
    G : dimension 1
    Mu : N_comps-1
    Hessian : (N_comps-1)X(N_comps-1) 
    
    For each phase for each of Mu,G,Hessian one .c is written
    
    Input : 
        
        phase_equil : phases for which derivative is evaluated     
    
    return :
        
        no return
    
    """
    
    print('Generating .c files')       
    print('')
    
    for phase in phase_equil:        
    
        phase_mod = Model(db_sys, component_sys, phase)
    
        my_list = [v.Y(phase,0,comp_calc[i]) for i in range(len(comp_calc))]
        my_list.append(v.T)
        print(my_list)
        
        GM_phase = phase_mod.GM
        
        y = sym.MatrixSymbol('y', phase_mod.GM.args.__len__()-1, 1)
        T = sym.MatrixSymbol('T', 1, 1)
        var = tuple(y)+tuple(T)
        print(var)
        
        mapped = zip(my_list, list(var))
        mapped = list(mapped)
        state_array_map= dict(zip(my_list, list(var)))
        
        GM_phase = phase_mod.GM.xreplace(state_array_map)    
        Mu_phase = {}
        Hess_phase = {}
        
        for n in var[0:-2]:
       
            Mu_phase[n] = GM_phase.diff(n) - GM_phase.diff(var[-2])
           # print(Mu_phase[n])
            
        for i in var[0:-2]:
            Hess_phase[i] = {}
            for j in var[0:-2]:
                Hess_phase[i][j] = Mu_phase[i].diff(j) - Mu_phase[i].diff(var[-2]) 
        
        G = sym.Symbol('G')
        
        Mu_phase = [Mu_phase[n] for n in var[0:-2]]
        Mu = sym.MatrixSymbol('Mu', N_comps-1,1)
        Mu_expr = sym.Matrix(Mu_phase)
        
        Hess_expr = [[Hess_phase[n1][n2] for n2 in var[0:-2]] for n1 in var[0:-2]] 
        dMudc = sym.MatrixSymbol('dMudc', N_comps-1, N_comps-1)
        dMudc_expr = sym.Matrix(Hess_expr)
           
        eqn_G = sym.Eq(G, GM_phase)
        eqn_Mu = sym.Eq(Mu, Mu_expr)
        eqn_Hess = sym.Eq(dMudc, dMudc_expr)
        
        if phase == Liquid_Phase:
       
           strt_G = 'G_LIQUID' 
           strt_Mu = 'Mu_LIQUID'  
           strt_Hess = 'Hess_LIQUID' 
              
        else:
          
           strt_G = 'G_SOLID' 
           strt_Mu = 'Mu_SOLID'  
           strt_Hess = 'Hess_SOLID' 

        codegen((strt_G, eqn_G), language='c', to_files=True)
        codegen((strt_Mu, eqn_Mu), language='c', to_files=True)
        codegen((strt_Hess, eqn_Hess), language='c', to_files=True)

    return
    
# In[]:        
#gen_c_phase_eq(phases)  
  
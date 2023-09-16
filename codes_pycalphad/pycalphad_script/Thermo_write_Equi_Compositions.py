#!/usr/bin/env python
# coding: utf-8

# In[]:
 
import sys    
 
import numpy as np

from pandas import DataFrame 

from pycalphad import Model, Database, equilibrium
import pycalphad.variables as v

import sympy as sym

from sympy.utilities.codegen import codegen

import Liquidus_Temp as Solv_Liq_Temp 

# In[]:

"""
 Script to write equilibrium compositions and second derivatives of free
 Energy to a .txt file and to generate .c files of free energy and derivatives

 This .py Script uses Temp_Range_Equilibrium() from Liquidus_Temp_Final
 imported as Solv_Liq_Temp   
 
 Input for Temp_Range_Equilibrium: 
    
    Alloy Composition, phases and Solid phase in equilibrium with Liquid Phase

 Input from User for the script:
     
          1: Database file(.tdb) 
          2: Solid Phase  
          3: Liquid Phase
          4: Component for Calculation 
          5: Alloy Composition    
          
"""

# In[]: 

db_sys = Database(sys.argv[1])        # Input no. 1
# Ex:  db_sys = Database("alzn_mey.tdb")

phases_sys = list(db_sys.phases.keys())
phases_sys = sorted(phases_sys) 
print('phase Present in system are', phases_sys)
print('')

Solid_Phase = sys.argv[2]  # Input no. 2
Liquid_Phase = sys.argv[3] # Input no. 3

# Ex : Solid_Phase = "FCC_A1"
# Ex : Liquid_Phase = "LIQUID"

phases_calc = [Solid_Phase,Liquid_Phase]  # Phases for Calculation  
phases_calc_sort = sorted(phases_calc)    # Phases Sorted       

print('Solid Phase :', phases_calc[0])
print('Liquid Phase :', phases_calc[1])
print('')

component_sys = list(db_sys.elements)
for component in component_sys:

    if component == '/-':
        component_sys.remove(component)               

print('Components in system are', component_sys)
print('')

# In[]:

calc_component_sys = [x for x in component_sys if x != 'VA']
calc_component_sys.sort()

Compo_calc_1 = sys.argv[4]  # Component for Calculation  Input no. 4

# Ex : Compo_calc_1 = 'ZN'  # Component for Calculation

ind_ex_1 = calc_component_sys.index(Compo_calc_1)    

print('Component for calculation :' , Compo_calc_1)
print('')

index_arr = np.arange(0,len(calc_component_sys),1)

ind_ex_2 = [x for j,x in enumerate(index_arr) if j!=ind_ex_1][0]  
Compo_calc_2 = calc_component_sys[ind_ex_2]

print('Other Component :' , Compo_calc_2)
print('')

# In[]:

c_alloy = float(sys.argv[5])   # Input no. 5

#c_alloy = 0.4
Temp_Tol = 0.001
Comp_Tol = 0.0005

print('Alloy Composition :', c_alloy)
print('')
# In[]:

strt_A = 'X_' + Compo_calc_1 + '_' + phases_calc[0]
strt_B = 'X_' + Compo_calc_1 + '_' + phases_calc[1]
strt_C = 'd2G_' + Solid_Phase + '_' + Compo_calc_1 + '_' + str(c_alloy)
strt_D = 'd2G_' + Liquid_Phase + '_' + Compo_calc_1 + '_' + str(c_alloy)


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
    
    print(phase)
    
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

def d_dG_eq(phase):
    
    """ 
    
    Function that returns double derivative of free energy function with 
    respect to component for calculation (Compo_Calc_1) in symbolic form
    
    Input : 
        
        phase : phase for which derivative is evaluated 
    
    return : 
        
        Symbolic form of double derivative
    
    """
    
    phase_mod = Model(db_sys, component_sys, phase)
    var = [v.Y(phase,0,calc_component_sys[ind_ex_1]), v.Y(phase,0,calc_component_sys[ind_ex_2]), v.T]
          
    GM_phase = phase_mod.GM
    GM_phase_diff = {}
      
    for n in var[0:-1]:
       
        GM_phase_diff[n] = GM_phase.diff(n)
        
    dG_phase = GM_phase_diff[var[0]] - GM_phase_diff[var[1]]
   
    Dou_dg_phase = dG_phase.diff(var[0]) - dG_phase.diff(var[1]) 

    return Dou_dg_phase 

# In[]:

#print(d_dG_eq(Solid_Phase)) 

# In[]:

def Tem_equi_composition(start,stop,c_alloy):
    
    """
    Function to write equi.compo(w.r.t Compo_Calc_1) and double derivatives 
    of each phase evaluated at equi. compositions at the Temperature 
    
    Inputs :
        
      start : temperature at which equilibria starts for the alloy composition 
      stop : temperature at which equilibria ends for the alloy composition   
      c_alloy: alloy composition: c_alloy
      
    return:   
      
      dataframes Containing outputs in order shown below 
      
      Dataframe 1 : Temp , equi. composition Solid, equi. composition Liquid   
      Dataframe 2 : Temp , double derivative 
      Dataframe 3 : Temp , equi. composition Solid, equi. composition Liquid
      
    """
    
    del_ta = 1
    
    df_Final = DataFrame()
    df_der_Solid = DataFrame()
    df_der_Liquid = DataFrame()

    T_array = np.arange( start, stop+ del_ta, del_ta) 
    
    for Temp in T_array:        
        
        phase_composition_equi = []
        
        eq1 = equilibrium(db_sys, component_sys, phases_sys, {v.X(Compo_calc_1): c_alloy ,v.T: Temp, v.P: 101325},output='GM')
        Phase_equili_T = np.array(eq1.Phase.squeeze())
        Phase_equili_T = [x for x in Phase_equili_T if x]
        Phase_equili_T = sorted(Phase_equili_T)
        
        df_Temp = DataFrame()
        der_1 = DataFrame()   
        der_2 = DataFrame()
    
        if phases_calc_sort == Phase_equili_T :                        
            
            for phase_equi in phases_calc:    
            
            #phases_calc : Ist element is always Solid and last is Liquid
              
                phase_compo_X = np.array(eq1.X.where(eq1.Phase == phase_equi).sel(component = Compo_calc_1).squeeze())
                phase_compo_X = phase_compo_X[~(np.isnan(phase_compo_X))] 
                phase_composition_equi.append(phase_compo_X)
                
                if phase_equi == Solid_Phase:
                   
                   my_list = [v.Y(phase_equi,0,calc_component_sys[ind_ex_1]), v.Y(phase_equi,0,calc_component_sys[ind_ex_2]), v.T]
                   der_val = d_dG_Sol.subs([(my_list[0], phase_compo_X[0]),(my_list[1], 1-phase_compo_X[0]),(my_list[2], Temp)])
                   der_1 = DataFrame({'T':Temp, strt_C: der_val },index=[0])                
                
                if phase_equi == Liquid_Phase:
                
                   my_list = [v.Y(phase_equi,0,calc_component_sys[ind_ex_1]), v.Y(phase_equi,0,calc_component_sys[ind_ex_2]), v.T]
                   der_val = d_dG_Liq.subs([(my_list[0], phase_compo_X[0]),(my_list[1], 1-phase_compo_X[0]),(my_list[2], Temp)])
                   der_2 = DataFrame({'T':Temp, strt_D: der_val },index=[0])           
   
                
            df_Temp = DataFrame({'T':Temp, strt_A: phase_composition_equi[0],strt_B: phase_composition_equi[1]},index=[0])
            df_Final = df_Final.append(df_Temp)
                 
            df_der_Solid = df_der_Solid.append(der_1)
            df_der_Liquid = df_der_Liquid.append(der_2)
                            
    return df_Final, df_der_Solid, df_der_Liquid



# In[]:

def gen_c_phase_eq(phase_equil):
     
    """
    Function to generate c codes with Free Energy(G), derivative of G 
    and double derivative w.r.t compo_calc_1 for each phase in phase_equil     
    
    Input : 
        
        phase_equil : phases for which derivative is evaluated     
    
    return :
        
        no return
    
    """
    
    print('Generating .c files')       
    print('')
    
    for phase in phase_equil:        
    
        phase_mod = Model(db_sys, component_sys, phase)
                        
        y = sym.MatrixSymbol('y', phase_mod.GM.args.__len__()-1, 1)
        T = sym.MatrixSymbol('T', 1, 1)
        var = tuple(y)+tuple(T)
        my_list = [v.Y(phase,0,calc_component_sys[ind_ex_1]), v.Y(phase,0,calc_component_sys[ind_ex_2]), 'T']
               
        mapped = zip(my_list, list(var))
        mapped = list(mapped)
        state_array_map= dict(zip(my_list, list(var)))
            
        GM_phase = phase_mod.GM.xreplace(state_array_map)
        
        GM_phase_diff = {}
           
        for n in var[0:-1]:
            
            GM_phase_diff[n] = GM_phase.diff(n)
        
        dG_phase = GM_phase_diff[var[0]] - GM_phase_diff[var[1]]
        
        Dou_dg_phase = dG_phase.diff(var[0]) - dG_phase.diff(var[1]) 
        
        G = sym.Symbol('dG')
        dG = sym.Symbol('dG')
        d2_G = sym.Symbol('d2_G')
        eqn_G = sym.Eq(G, GM_phase)
        eqn = sym.Eq(dG, dG_phase)
        eqn1 = sym.Eq(d2_G, Dou_dg_phase)
        
        if phase == Liquid_Phase:
          
            strt_ther = 'G_' +  phase 
            strt_thermo = 'DG_' +  phase  
            strt_chempot = 'dou_DG_' + phase
    
        else:
          
            strt_ther = 'G_SOLID' 
            strt_thermo = 'DG_SOLID'  
            strt_chempot = 'dou_DG_SOLID' 
    
        codegen((strt_ther, eqn_G), language='c', to_files=True)
        
        codegen((strt_thermo, eqn), language='c', to_files=True)

        codegen((strt_chempot, eqn1), language='c', to_files=True)


    return
# In[]:

d_dG_Sol = d_dG_eq(Solid_Phase)
d_dG_Liq = d_dG_eq(Liquid_Phase)
    
T_Range = Solv_Liq_Temp.Temp_Range_Equilibrium(c_alloy,phases_sys,Solid_Phase)

T_Range = sorted(T_Range)
start = T_Range[0]

if (start - int(start)) != 0:
    
    start = int(start) + 1

stop = int(T_Range[1])

df_compo, d2_G_Sol, d2_G_Liq = Tem_equi_composition(start,stop,c_alloy)

# In[]:

df_compo.to_csv("Equilibrium_Composition.txt", sep = ' ',index = False)
d2_G_Sol.to_csv(strt_C, sep = ' ',index = False)
d2_G_Liq.to_csv(strt_D, sep = ' ',index = False)

# In[]:
    
gen_c_phase_eq(phases_calc)

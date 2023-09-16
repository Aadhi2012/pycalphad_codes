#!/usr/bin/env python
# coding: 
# In[]:

import sys

import numpy as np

from pycalphad import Database, equilibrium
import pycalphad.variables as v


# In[]:

"""
  
 This .py script contains functions for the following purposes:
     
 Temp_Range_Equilibrium:
     
     To find the range of temperature over which an equilibrium exists for 
     a given two phase : Here Solid and Liquid Phase for an alloy composition
     
 Solid_Liq_Equilibrium :     

     For a given alloy composition, the solid phase in equilibrium with the
     liquid phase is found and the corresponding temperature range over which 
     equilibrium exists
     
"""
# In[]: 

db_sys = Database(sys.argv[1])        # Input no. 1

#db_sys = Database("alzn_mey.tdb")

phases_sys = list(db_sys.phases.keys())
phases_sys = sorted(phases_sys) 

Solid_Phase = sys.argv[2]  # Input no. 2
Liquid_Phase = sys.argv[3] # Input no. 3

#Solid_Phase = "FCC_A1"
#Liquid_Phase = "LIQUID"

component_sys = list(db_sys.elements)
for component in component_sys:

    if component == '/-':
        component_sys.remove(component)               

#print('Components in system are', component_sys)
print('')

# In[]:

calc_component_sys = [x for x in component_sys if x != 'VA']
calc_component_sys.sort()

Compo_calc_1 = sys.argv[4]  # Component for Calculation  Input no. 4

#Compo_calc_1 = calc_component_sys[1]  # Component for Calculation
ind_ex_1 = calc_component_sys.index(Compo_calc_1)    


# In[]:

Temp_Tol = 0.0001      # Temperature tolerance for searches
Comp_Tol = 0.00005     # composition tolerance
T_step = 0.0001        # For finding slope 
    
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

# In[]:

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

def Binary_Search(T_Low,T_High,c_alloy,phases_in_range,phases_calc_sort): 
       
    """
    
    Binary Search to find One Temp where Liquid,Solid equilibrium exists 
    
    Inputs:
        
        T_Low : Temperature above which Solid exists                                                     
        T_High : Temperature above which Liq exists 
    
    Between T_High and T_Low , (Liq + Solid) exists found using Binary Search 
    
    
        phases_in_range :phases considered to find temperature range 
                         of equilibrium
    
       
        phases_calc_sort : phases for calculation [Solid,Liquid] sorted
            
    Search stops when Temp_diff acheives tolerance or Temp for equilibrium
    is found(flag = 1)  
    
    return:
        
        Function returns an array containing following entries:
        
        Array : [T_Low, T_current, T_High, flag]
    """
    
    flag = 0
    Temp_diff = abs(T_High - T_Low)
    
    T_Current = (T_High + T_Low)*0.5
    phase_equi_cur = phase_in_equilibrium(T_Current,c_alloy,phases_in_range)      
    
    if phase_equi_cur == phases_calc_sort:
        
        flag = 1
                                                                                                                
    while(flag == 0 and Temp_diff >= Temp_Tol):  
        
        phase_equi_cur=phase_in_equilibrium(T_Current,c_alloy,phases_in_range)
        
        if phase_equi_cur == phases_calc_sort:
        
            flag = 1
            break
        
        elif phase_equi_cur == [Liquid_Phase]:
                
            T_High = T_Current
        
        else:
        
            T_Low = T_Current
        
        T_Current = (T_High + T_Low)*0.5
        Temp_diff = abs(T_High - T_Low)     
    
    if(flag == 0 and Temp_diff < Temp_Tol):
        
        print('Temperature for Equilibrium not found  :' ,phases_calc_sort )
        
        
    return [T_Low,T_Current,T_High,flag]                

# In[]:

def Binary_Melting(T_Low,T_High,c_alloy,phase_tem,phases_in_range,phases_calc_sort):
    
    """
    Functions to find Liquidus Temp (phase_tem is Liquid)
    Functions to find Solidus Temp (phase_tem is Solid)
    
    Inputs:
    
    T_High : phase Tem : Solid,  Temperature at which (Liq + solid) exists 
    T_Low : phase Tem : Solid,  Temperature at which solid exists 
    
    T_High : phase Tem : Liquid,  Temperature at which Liquid exists 
    T_Low : phase Tem : Liquid,  Temperature at which (Liq + solid) exists 
    
    phases_in_range : phases considered to find temperature range
    
    Between T_High and T_Low , (Liq + Solid) exists found using Binary Search 
    till Temp_tol is acheived. Search stop when Temp_diff acheives a tolerance 
    and Equilibria exists at the Temp  found(flag = 1)  
    
    return:
        
    Solidus or Liquidus Temperature    
    
    """
    
    flag = 0
    
    if phase_tem == [Liquid_Phase]:
       
        T_Current = T_Low
    
    else:
          
        T_Current = T_High  
        
    phase_equi_cur = phase_in_equilibrium(T_Current,c_alloy,phases_in_range)
    
    if phase_equi_cur == phases_calc_sort:
        
        flag = 1
    
    Temp_diff = abs(T_High - T_Low)
    
    while(not(bool(abs(Temp_diff) < Temp_Tol and flag == 1))):  
        
        T_Current = (T_High + T_Low)*0.5
        phase_equi_cur=phase_in_equilibrium(T_Current,c_alloy,phases_in_range)
        
        if phase_equi_cur == phases_calc_sort:
            
            flag = 1
            
            if phase_tem == [Liquid_Phase] :
                
                T_Low = T_Current
                
            else:
                
                T_High = T_Current
            
        else:
           
            flag = 0
            
            if phase_tem == [Liquid_Phase] :
                
                T_High = T_Current
                
            else:
                
                T_Low = T_Current
                
        Temp_diff = abs(T_High - T_Low)
      
    if(flag == 0 and Temp_diff < Temp_Tol):
        
        print('Temperature not found')
                
    return T_Current                


# In[]:

def Temp_Range_Equilibrium(c_alloy,phases_in_range,Solid_Phase):
    
    """
    Function to find Temp range over which Liq + Solid equilibria exists
    Liquidus and Solidus Temp for an alloy composition is found
    
    Inputs:
        
    c_alloy: Alloy Composition    
    phases_in_range : phases considered to find temperature range
    Solid_Phase : Solid phase in equilibrium with Liquid
    
    return:
        
        Returns range of temperature for equilibrium
        
        Output : [Solidus Temperature, Liquidus Temperature]
    
    """
    
    start = 400       # start : Temperature below which only Solid Exists
    stop = 1800       # stop : Temperature above which only Liquid Exists
    
    phases_calc_sort = sorted([Solid_Phase,Liquid_Phase])
    
    T_Bin_Res = Binary_Search(start,stop,c_alloy,phases_in_range,phases_calc_sort) 
    T_Bin_Res[0:-1] = sorted(T_Bin_Res[0:-1])

    if T_Bin_Res[3] == 1: 
        
        print('Found_Equilibrium : ',phases_calc_sort)
        print('')
        T_Liquidus=Binary_Melting(T_Bin_Res[1],T_Bin_Res[2],c_alloy,[Liquid_Phase],phases_in_range,phases_calc_sort)
        
        print('Liquidus Temperature : ',T_Liquidus)
        print('Liq Composition : ',equi_composition_phase(T_Liquidus,c_alloy,Liquid_Phase,phases_in_range))
        print('')
        T_Solidus=Binary_Melting(T_Bin_Res[0],T_Bin_Res[1],c_alloy,[Solid_Phase],phases_in_range,phases_calc_sort)    
        
        print('Solidus Temperture : ',T_Solidus)
        print('Solid Composition : ',equi_composition_phase(T_Solidus,c_alloy,Solid_Phase,phases_in_range))
        print('') 
        return [T_Solidus,T_Liquidus]

    else:

        print('Equilibrium Not Found : ',phases_calc_sort)


        return [0,0]

# In[]:

def Solid_Liq_Equilibrium(C_0, Solid_Phases):
    
    """
    Function to find temp range of two phase equilibria and the phases in 
    equilibrium
    
    Inputs: 
        
          Alloy Composition : C_0 w.r.t.component for calculation
          Solid_Phases : Array containing Solid phases in the system 
    
    returns:
          
            Function returns array Temperature range and phases in equilibrium 
            [[Solidus Temperature, Liquidus Temperature],[Solid, Liquid]] 
            respectively for each equilibria present
            
    """
    res_arr = []
    
    for Solid_Phase in Solid_Phases:

        T_res = Temp_Range_Equilibrium(C_0,phases_sys,Solid_Phase)

        if T_res == [0,0]:

            continue

        else:
           
            Solid_in_equi = Solid_Phase
            phases_equi = [Solid_in_equi,Liquid_Phase]
            
            app_arr = np.array([T_res,phases_equi])
            res_arr.append(app_arr)
    
    res_arr = np.array(res_arr)
    
    return res_arr  


# In[]:

def sign_slope_Liquidus(C,T,phases_calc_sort):
    
    """
    Function to find Liquidus slope
    
    Inputs:
    
        C : Composition to find Liquidus line slope at Temperature T 
        At ( C,T ) Solid + Liquid equilibria should exist
        phases_calc_sort : phases in Equilibrium ( Solid + Liquid )
    
        c_perturb: A point on Liquidus line has equilibrium phase LIQUID,hence 
                   to find Composition of solid c_perturb is found
    
    """
    
    global c_perturb
    
    T_for = T - T_step
    c_calc = C
    phases_C = phase_in_equilibrium(T,c_calc,phases_sys)

    if phases_C == phases_calc_sort:
    
        Liq_compo_1 = equi_composition_phase(T,c_calc,Liquid_Phase,phases_sys)[0]
        Liq_compo_2 = equi_composition_phase(T_for,c_calc,Liquid_Phase,phases_sys)[0]    
    
        Liq_slope = (Liq_compo_1 - Liq_compo_2)/T_step
        
        if Liq_slope > 0:

            c_perturb = 0.0005

        else:

            c_perturb = -0.0005    
            
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
    
    Liq_compo_1=equi_composition_phase(T,c_slope+c_perturb,Liquid_Phase,phases_sys)[0]
    Liq_compo_2=equi_composition_phase(T_for,c_slope+c_perturb,Liquid_Phase,phases_sys)[0]    
    
    slope_Liq = T_step/(Liq_compo_1 - Liq_compo_2)
    
    x = new_x(slope_Liq, T, Liq_compo_1 - c_alloy)
    
  #  print(x)
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
    
    c_alloy : Alloy Composition
    c_sol_eval : half of Composition of solid + liquid
    C_eval: composition for which Liquidus temperature to be found
    T_new: Temperature at which the composition for calculation has to be found
    phases_calc_sort:  [Solid,Liquid] in equilibrum sorted
    
    """
    
    ph_alloy = phase_in_equilibrium(T_new,c_alloy, phases_sys) 
    ph_eval = phase_in_equilibrium(T_new,C_eval, phases_sys)
    ph_sol = phase_in_equilibrium(T_new,c_sol_eval, phases_sys)
    
    Solid_Phase =  [x for j,x in enumerate(phases_calc_sort) if j!= Liquid_Phase][0]
    
    ph_all = [ph_alloy,ph_eval,ph_sol]
    comp_all = [c_alloy,C_eval,c_sol_eval]
    
    C_S = -1
    C_L = -1
    
    for i in range(len(ph_all)):
        
        if ph_all[i] == phases_calc_sort:
              
            c = comp_all[i]
            return c
    
    for i in range(len(ph_all)):
        
        if ph_all[i] == [Solid_Phase]:
    
            C_S = comp_all[i]
        
        if ph_all[i] == [Liquid_Phase]:
    
            C_L = comp_all[i]
     
    if (C_S != -1 and C_L != -1):
        
        c = Compo_Binary_Search(C_S,C_L,T_new,phases_calc_sort)
        return c
        
    else:
        
        print('Composition for equilibrium not found')
        
        
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
        
        if (abs(C_S - C_L) < Comp_Tol):
            
            print('Composition for Equilibrium not found')
            break
        
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
    Solid_Phase =  [x for j,x in enumerate(phases_calc_sort) if j!= Liquid_Phase][0]
     
    C_eval = c_known
    T_eval = T_known
    c_sol_eval = 0.5*(equi_composition_phase(T_known,C_eval+ c_perturb, Solid_Phase, phases_sys)[0] + C_eval)
    
    comp_diff = c_alloy - c_known
    
    if abs(comp_diff) <= Comp_Tol:
                
        T_new = T_eval 
            
    while(abs(comp_diff) > Comp_Tol) :   
        
        T_new = x_Secant(c_alloy,C_eval,T_eval) 
 #       print(T_new)
        C_eql = find_c_equilibrium(c_alloy,c_known,c_sol_eval,T_new,phases_calc_sort)
        
        if C_eql != -1:
        
            c_liq_new = equi_composition_phase(T_new,C_eql,Liquid_Phase, phases_sys)
            
            comp_diff = c_liq_new - c_alloy
            T_eval = T_new
            C_eval = c_liq_new
        
        else:
            
            break
        
    phase_equi_cur = phase_in_equilibrium(T_new,c_alloy+ c_perturb, phases_sys)
    
    if phase_equi_cur == phases_calc_sort:
            
        #print('Liquidus Temperature Found : ' ,T_new)     
        return T_new,phase_equi_cur
        
    else:
           
        print('Liquidus Temperature not Found')             
        
        return T_new,phase_equi_cur 

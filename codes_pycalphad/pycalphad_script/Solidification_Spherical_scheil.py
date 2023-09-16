#!/usr/bin/env python
# coding: utf-8

"""

Script to simulate solidification in a Liquid droplet. 
Heat balance considering heat of fusion and convective transfer of heat
Log Normal Nucleation probability
Growth rate : Dendrite tip velocity

To obtain thermal history of the droplet for binary alloy system 

"""

# In[]:

from pycalphad import Database, equilibrium
import pycalphad.variables as v
from pandas import DataFrame as df

import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter 

import scipy.special as sp
import scipy.optimize as opti
from scipy.interpolate import interp1d, interp2d
import re    

# In[]

import Gibb_Thomson_Coef as GTC
import Liquidus_Tem as Solv_Liq_Temp          # Melting Point of Alloy
import Liquidus_Tem_Entire as Liq_Tem_ent     # Liquidus line

# In[]:

"""    

Inputs read from .in input file taken as argument using sys.argv

tdbfname, NUMCOMPONENTS, COMPONENTS, tdb_phases, num_thermo_phases, ceq, V 

"""
print('Reading Input file')

#with open(sys.argv[1], 'r') as my_file:    
with open('Input.in', 'r') as my_file:    
    
    flines=my_file.readlines()

    flg_tdbf, flg_Comps, flg_Num_Comps, flg_Phase, flg_Num_Phases = 0,0,0,0,0  
    flg_V = 0
    
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
                
        searchstatus = re.match(r'\bV\b', line)
        if searchstatus:
            flg_V = 1
            Vline = line
            V_m = Vline.split("=")[1].replace(';', '').strip()
            V_m = float(V_m)

        
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
    
#print(equi_composition_phase(793, Liq_Comp, Liquid_Phase, phases_sys))

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

print(' Alloy Composition :',Liq_Comp)

T_melt, Output_arr = Solv_Liq_Temp.Temp_Range_Equilibrium(350, 1800, Liq_Comp, phases_sys)

print(T_melt)

# In[]:

# slope for finding c_perturb for pyclaphad equilibrium calculations    
    
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
T_Liq = Liquid_file['Tem']
str_name = 'X_'+ Comp_Indep[0] # only for binary
Compo_Liq_arr = Liquid_file[str_name]

LiqT_arr = [ [Compo_Liq_arr[i], T_Liq[i]] for i in range(len(T_Liq)) ]  
LiqT_arr = np.array(Liq_T_arr) 
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

# In[]:

# c_liq_max and c_liq_min are found for the interpolation functions.
# The interpolation is done for liq composition. The corresponding solid 
# composition for the min or max of Liq Composition will be outside the 
# interpolation range for finding melting point      

c_max = np.max(Compo_Liq_arr) - del_c
c_min = np.min(Compo_Liq_arr) + del_c

c_liq_max, c_liq_min = c_max, c_min

c_sol_min = f_cubic_Liq_Sol(c_min) 
c_sol_max = f_cubic_Liq_Sol(c_max) 

if c_sol_max > c_max :
    
    c_liq_max = f_cubic_Sol_Liq(c_max)
    
elif c_sol_min < c_min:
      
    c_liq_min = f_cubic_Sol_Liq(c_min)
       
print(c_liq_min,c_liq_max)     

# In[]:
     
Ga_mma = GTC.Gibbs_Thom_Coef(350, 1800, Liq_Comp) #m-K Gibbs Thomson Coefficient
D = 1*10**-9       #m**2/s  # Diffusivity
error_tolerance = 5e-4   # for fsolve errors 

# In[62]:

def solve_eqn(G,C_0,T):
    
    # Equation from LGK formulation for obtaining V tip and R tip
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
   # print(c_l,c_s)
    Tm_cs = f_cubic_Liquidus(c_s).item()
    
    d_tip = abs(Ga_mma/(Tm_cs - Tm_cl))*x_tip
    eq_1 = abs(2*Ga_mma*x_tip - (Tm_cl - T))
    eq_2 = abs(Pe - 2*(np.pi)**2*d_tip)
    eq_3 = abs((c_l - C_0)/(c_l - c_s) - Pe*(np.exp(Pe))*sp.exp1(Pe))
   
    return[eq_1,eq_2,eq_3] 
 
# In[]:    
   
def find_V_R_tip(Guess,C_far,T):
    
    # Find V and R of dendritic tip using solv_eqn
    
    res_T, info, ier, msg = opti.fsolve(solve_eqn,Guess,args = (C_far,T),full_output= True)
    res_T = abs(res_T)   # absolute value is the solution
    Guess = res_T      # for next iteration updated Guess is used
        
    error = solve_eqn(res_T,C_far,T) 
    error = np.array(error)

    R_tip = 1/res_T[1]    
    V_tip = 2*D*res_T[2]*res_T[1]
    
    if error.max() > error_tolerance:

        print('error greater than error_tolerance :', error_tolerance)
        
        return V_tip, R_tip, res_T
    
    #print( ' error :' , error)
     
    return V_tip, R_tip, res_T

    
# In[]:

def find_ME_phase(compo,Tem, Phase):
    
    """
    
     Find Molar enthalpy of phase : unit J/mol
     Molar Enthalpy for a phase at Temperature is found using equilibrium
     
     Inputs:
         Compo : Composition
         T : Temperature
         Phase : Phase for which Molar enthalpy is found
         
     Output:
         
         Molar Enthalpy in units J/mol
     
    """
    
    conds = {v.T: Tem, v.P:101325}
    var_comp = [v.X(comp_calc[i]) for i in range(len(comp_calc)-1)]
    conds_comp = {var_comp[i]: compo[i] for i in range(len(var_comp))}
    
    conds.update(conds_comp)
    
    cal_M = equilibrium(db_sys, component_sys, [Phase], conds,output='HM')
    ME_Phase = cal_M.HM.squeeze()
    
    return ME_Phase 

# In[]:
    
#print(find_ME(1 - 1e-7,694,Liquid_Phase) - find_ME(1 - 1e-7,694,'HCP_A3'))

# In[]:

def find_C_phase(compo,Tem, Phase):
    
    """
    Find Heat capacity for a phase at a temperature
    
    Inputs:
         Compo : Composition
         T : Temperature
         Phase: Phase for which heat capacity is calculated
         
     Output:
         
         Heat Capcity in units J/mol-K
     
    """
    
    conds = {v.T: Tem, v.P:101325}
    var_comp = [v.X(comp_calc[i]) for i in range(len(comp_calc)-1)]
    conds_comp = {var_comp[i]: compo[i] for i in range(len(var_comp))}
    
    conds.update(conds_comp)
    
    cal = equilibrium(db_sys, component_sys, [Phase], conds,output='heat_capacity') 
    heat_cap = cal.heat_capacity.squeeze().item()
    
    return heat_cap 

# In[]:

# Heat of fusion can be taken as a constant calculated at the melting point 
# of the alloy composition 
    
C_Liq_eq = Liq_Comp
C_Sol_eq = equi_composition_phase(T_melt, Liq_Comp, Solid_Phase, phases_sys)[:,0]

H_L = find_ME_phase(C_Liq_eq, T_melt, Liquid_Phase)
H_S = find_ME_phase(C_Sol_eq, T_melt, Solid_Phase)

L_f = np.array(H_L - H_S)

# In[]:

# Heat of fusion can also be interpolated as a function of Temperature and 
# Composition as done below

"""     
d_com = 0.01  
c_mesh_arr = np.arange(d_com,1,d_com)
T_mesh_arr = np.arange(T_melt-100,T_melt+5,5)     

ME_S_arr = np.zeros((len(T_mesh_arr),len(c_mesh_arr)))
ME_L_arr = np.zeros((len(T_mesh_arr),len(c_mesh_arr)))

i  = 0 

for TT in T_mesh_arr:
    j = 0
    for cc in c_mesh_arr:
        
        ME_L_arr[i,j] = find_ME_phase([cc], TT, Liquid_Phase) 
        ME_S_arr[i,j] = find_ME_phase([cc], TT, Solid_Phase) 
        j += 1
        
    i +=1    

f_ME_L = interp2d(c_mesh_arr, T_mesh_arr, ME_L_arr, kind = 'cubic')
f_ME_S = interp2d(c_mesh_arr, T_mesh_arr, ME_S_arr, kind = 'cubic')
"""

# In[]:
    
def find_prob_nuclei(T,c_liq,sig):
    
    """
     Calculate probability of nucleation.
     Probability of nucleation given by log normal distribution  
    
    T : Temperature to calculate undercooling
    c_liq : Liquid alloy composition 
    sig : Standard deviation of nucleation probability
    """
    
    Tm_Cl = f_cubic_Liquidus(c_liq).item()
    T_cool = Tm_Cl - T 
    dT_m = 1
    dT_sig = sig
    
    pro_bab = 1/(2*np.pi)**0.5/dT_sig/T_cool*np.exp(-0.5*((np.log(T_cool)- np.log(dT_m))/dT_sig)**2)
 
    return pro_bab

# In[]:

# plot nucleation probability distribution

sig_arr = [0.5,0.75,1]    
    
for sig in sig_arr:    

    u_c_arr = np.arange(0.01,5.1,0.01)
    pr_arr = np.zeros(len(u_c_arr))
        
    for i in range(len(u_c_arr)): 
    
        pr_arr[i] = find_prob_nuclei(T_melt - u_c_arr[i], Liq_Comp, sig)
    
    plt.plot(u_c_arr,pr_arr, linewidth= 4, label = '[$\sigma $] ='+str(sig))
    plt.title('Nucleation Probability')
    plt.ylabel('probability')
    plt.xlabel('Undercooling')
    plt.legend()

# In[]:

def active_nucleation(T,c_liq, fract_solid, sig_ma):
  
    """  
    Check for occurrence of nucleation. 
    Using log normal distribution for nucleation probability
    
    Generate random number between 0 and 1 
    if random no less than prob then nucleate 
    
    """
    
    prob = find_prob_nuclei(T,c_liq, sig_ma)
    prob = prob*(1- fract_solid)
    n_s = 0
    rand_pro = np.random.random()
    if rand_pro < prob:
        n_s = 1
     
    return n_s

# In[]:

def frac_sol(r_array):

    # r_array : array containing radius of nuclei present at the time step
    # return fraction of solid    

    rad_cube = 0
    
    for r in r_array:    
        rad_cube += r**3        
    
    frac_sol = rad_cube/R_d**3  
    
    return frac_sol

# In[]:

def create_r_array(r_entire_array):    

    # r_entire_array contain radius with time for each nuclei present    
    # r_array contain radius of nuclei present at the time step
    
    le_n_entire = len(r_entire_array)
    r_array = np.zeros(le_n_entire)
    
    for i in range(le_n_entire):
        
        temp_arr = r_entire_array[i] # 0 corresponds to 1st nuclei formed
        r_array[i] = temp_arr[-1]   # Last element is radius at time step
    
    return r_array

# In[]:
    
Lf_arr = [] # array to store heat of fusion

# In[]:    

def change_temp(r_array,T_dt,c_dt,c_liq_s,c_sol,v_tip):    
    
    # Temperature change according to heat balance from heat of fusion from 
    # solidification and convective transfer of heat via surface of droplet
     
    # Temperature dependent Heat of fusion
    #L = f_ME_L(c_liq_s, T_dt).item() - f_ME_S(c_sol, T_dt).item()   
    #Lf_arr.append(L)
     
    # Heat of fusion : Constant 
    L = L_f
    
    dT_net = -h*3*R_d**2*V_m*(T_dt - 300)  # due to heat from convection
    
    rad_cube = R_d**3
    sol_vol = 0
    
    for r in r_array:    
       
        sol_vol += r**2*v_tip  # solidified volume : multiply sol_vol with dt 
    
    dT_net += sol_vol*L*3      # heat change due to solidification
    
    # Entire droplet heated containing solid and Liquid  
    dT_net = dT_net*dt/h_c/(rad_cube) 

    return dT_net
    
# In[]:
    
def change_composition(r_array,c_dt,v_tip,c_sol):    
    
    # Liq composition change according to scheil
    # c_dt is the liq composition at the time step
    # c_sol composition of solid forming from c_dt
    # v_tip is rate of solidification
    
    rad_cube = 0
    sol_vol = 0
      
    for r in r_array:    
       
        rad_cube += r**3        
        sol_vol += r**2 
    
    sol_vol = sol_vol*3*v_tip*dt
    rad_cube = R_d**3 - rad_cube  # Comp change only for Liquid unlike Temp change
    
    d_c_dis = sol_vol*(c_dt - c_sol)/(rad_cube)
 
    c_dt = c_dt + d_c_dis
    
    return c_dt

# In[]:

def check_in_solid(pos_it, pos_arr, r_array):
    
    # check if pos_it : [x,y,z] is inside solid
    # pos_arr : position of nucleation of previous nucleations
    
    flag = 0 
    len_r_arr = len(r_array) 
    
    if len_r_arr == 0: # No solid yet
        
        flag = 1
    
    for i in range(len_r_arr):
        
        pos_n = pos_arr[i]
        rad = r_array[i]
        ch_val = (pos_it[0] - pos_n[0])**2 + (pos_it[1] - pos_n[1])**2 + (pos_it[2] - pos_n[2])**2 
        
        if ch_val < (rad/R_d)**2: 
            flag = 0  # inside atleat one solid 
            break
        
        else:

            flag = 1
            
    return flag

# In[]:

def update_Radis_arr(r_ent_array,T_dt,c_dt,v_tip,n,time_nucl, pos_arr):    
     
    # Check for nucleation and add new nuclei
    """
    r_ent_array contains radius with time for each nuclei present  
    T_dt : Temperature of the droplet at the time step
    c_dt : Liquid composition at the time step
    time_nucl : contains time nucleation of nucleations so far
    pos_arr : contains nucleation sites so far
    v_dt : increase radius by v_tip*dt for solid present 
    n: Check for nucleation every n timesteps
    
    """
    
    r_array = create_r_array(r_ent_array)
    fra_sol = frac_sol(r_array) 
    
    Tm_Com = f_cubic_Liquidus(c_dt).item()
    u_c_dt = Tm_Com - T_dt  # Undercooling w.r.t c_dt
    
    for i in range(len(r_ent_array)):    
        
        r_num_array = r_ent_array[i] 
        # d_r_inc is the new radius after solidification
        d_r_inc = r_num_array[-1] + v_tip*dt  
        r_num_array = np.concatenate((r_num_array, np.array([d_r_inc])))
        r_ent_array[i] = r_num_array
        
    if (n+1)%10 == 0 and u_c_dt > u_c_thr:
          
       # u_c_thr is undercooling threshold: Above u_c_thr check for nucleation  
       # nucleation at a site where there is no solid 
 
       # checking for nucleation 
       d_new_nuclei = active_nucleation(T_dt, c_dt,fra_sol, sig_ma)
           
       if d_new_nuclei == 1 :
         
           x = np.random.random()
           y = np.random.random()
           z = np.random.random()
           
           pos_add = [x,y,z]  
           check_sol = check_in_solid(pos_add,pos_arr, r_array)
               
           while(check_sol == 0):
               
               # find x,y,z for nucleation
               
               x = np.random.random()
               y = np.random.random()
               z = np.random.random()
               pos_add = [x,y,z]
               check_sol = check_in_solid(pos_add,pos_arr, r_array)
           
           if check_sol == 0:
                print('')
         
           time_nucl.append(n)
           r_add = [0]
           pos_arr.append(pos_add)
           r_ent_array.append(r_add)
            
    return r_ent_array, time_nucl, pos_arr

# In[]:

x = 20
R_d = x*10**-6       # Radius of Droplet : in m
sig_ma = 0.5

h_c = find_C_phase(Liq_Comp,T_melt, Liquid_Phase) # heat capcity of Liquid
#print(h_c)

h = int(10**4*R_d*h_c/(T_melt - 300)/3/V_m) # heat transfer coefficient
A = int(h*V_m*3*(T_melt-300)/R_d/h_c) # cooling rate
#print(A)

# In[]:

h_arr =[int(h),h*0.5,h*0.25]    
cool_arr = [A, int(A*0.5),int(A*0.25)]       
th_arr = [0,0.3,1]               

# In[]:

V_tip_ar_1_ent = {}
# dT_1_ent : change in temp at every timestep,T_inter_ent: Interface Temp
dT_1_ent, T_inter_ent = {}, {} 
u_c_chan_ent = {}    # Undercooling
c_far_ent = {}   # far field composition
Radis_ent = {}   
frac_sol_ent = {}
timer_ent = {}

pos_arr_ent = {} # Nucleation position array
time_nucl_ent = {} # Nucleation time array

# In[]:
    
i = -1    
    
for h in h_arr:
    
    i += 1
    print(h)
    time_calc = 0
    n_time_step = 0
    num_step = 75000         # total time steps
    u_c_thr = th_arr[i]       
    dt = 3e-7                 # time step

    pos_arr = []
    time_nucl = []
    
    C_dis = Liq_Comp
    T_1 = T_melt
    
    print('C_initial, Tm  :',C_dis,T_1)
    
    dT_1 = -(h*3*R_d**2*V_m)/h_c/(R_d**3)*(T_melt - 300)*dt
    
    V_tip_ar_1 = []
    dT_1_arr, T_inter_arr = [], []
    u_c_chan_arr = []   
    c_far_arr = []  
    Radis_arr = []
    frac_sol_arr = []
    timer_arr = []
    
    Guess = np.append(C_dis,0.1) # Guess for finding dendrite tip velocity
    Guess = np.append(Guess,0.1)
    C_new = C_dis
    print('Start')
    
    while(n_time_step < num_step):
        
        T_1 += dT_1                       # Initial Undercooling
        C_dis = C_new                     # Liquid Composition
        dT_1_arr.append(dT_1)
        timer_arr.append(time_calc)    

        # find V tip and R tip using find_V_R_tip
        V_tip_sol , R_tip_sol, res_T  =  find_V_R_tip(Guess, C_dis, T_1)
    
        T_int_sol = T_1          # InterfaceTemp
        c_liq_sol = res_T[0]     # Comp Liq solved from find_V_R_tip 
        
        V_tip_ar_1.append(V_tip_sol)
        T_inter_arr.append(T_int_sol)
        
        # u_c_cur : current undercooling w.r.t C_dis
        u_c_cur = f_cubic_Liquidus(C_dis).item() - T_1 
        u_c_chan_arr.append(u_c_cur)
    
        c_sol_solve = f_cubic_Liq_Sol(c_liq_sol) # solid comp in eqb with Liquid
        c_far_arr.append(C_dis)
        
        Guess = res_T # Guess for next iteration
        
        r_array = create_r_array(Radis_arr)
        
        # New Temp and Liq comp for next iteration
        dT_1 = change_temp(r_array, T_1, C_dis, c_liq_sol,c_sol_solve, V_tip_sol)   
        C_new = change_composition(r_array, C_dis, V_tip_sol, c_sol_solve) 
        
        Radis_arr, time_nucl, pos_arr = update_Radis_arr(Radis_arr, T_1, C_dis, V_tip_sol, n_time_step, time_nucl,pos_arr)
        
        if (n_time_step)%1000 == 0:
            
             print(dT_1, u_c_cur, n_time_step, len(Radis_arr))
             print('f_s :' , frac_sol(r_array))
             print('c_s:', c_sol_solve,'c_l :', c_liq_sol)
             print("")
        
        sol_fr = frac_sol(r_array)
        frac_sol_arr.append(sol_fr)
        time_calc += dt 
        n_time_step += 1
    
         
    print(C_dis,T_1,u_c_chan_arr[-1], len(Radis_arr))    
    
    pos_arr_ent[h] = pos_arr
    time_nucl_ent[h] = time_nucl
    V_tip_ar_1_ent[h] = V_tip_ar_1
    dT_1_ent[h], T_inter_ent[h] = dT_1_arr, T_inter_arr
    u_c_chan_ent[h] = u_c_chan_arr   
    c_far_ent[h] = c_far_arr  
    Radis_ent[h] = Radis_arr
    frac_sol_ent[h] = frac_sol_arr
    timer_ent[h] = timer_arr
    
    print('done')    

# In[]:
    
print((time_nucl_ent[h_arr[0]]))    

# In[]:
    
print(len(Radis_arr[0]))

# In[]:
    
for i in range(len(h_arr)):    
    
    timer_plot_arr = np.zeros(len(time_nucl_ent[h_arr[i]]))    
    tem_int_plot_arr =  np.zeros(len(time_nucl_ent[h_arr[i]]))
    tem0p_arr =  time_nucl_ent[h_arr[i]]
    int_tem_arr = T_inter_ent[h_arr[i]]

    for j in range(len(timer_plot_arr)):
    
      timer_plot_arr[j] = dt*tem0p_arr[j]
      tem_int_plot_arr[j] = int_tem_arr[tem0p_arr[j]]

    plt.plot(timer_ent[h_arr[i]], T_inter_ent[h_arr[i]], label = ('h = '+str(h_arr[i])+' '+str(cool_arr[i])), linewidth = 3)
    plt.scatter(timer_plot_arr, tem_int_plot_arr, color='orange')
    plt.title('Tem vs time : ' + '$\sigma$  = '+ str(0.5))
    plt.ylabel('Temp')
    plt.xlabel('time')
    plt.legend()

# In[]:
    
#print(np.max(L_arr))    
    
# In[]:

for i in range(len(h_arr)):    
    
    plt.plot(timer_ent[h_arr[i]], dT_1_ent[h_arr[i]], label = ('h = '+str(h_arr[i])+' '+str(cool_arr[i])), linewidth = 3)
    plt.title('frac_sol vs Liq_comp : ' + '$\sigma$  = '+ str(0.5))
    plt.ylabel('Fraction of Solid')
    plt.xlabel('Liq_Comp')
    plt.legend()

plt.axhline(0, xmin = 0, color = 'red', linestyle = '--')  
  
# In[]:

#plot size distribution
    
cur_i = 2   

Radis_arr_cur = Radis_ent[h_arr[cur_i]]
time_nucl_cur = time_nucl_ent[h_arr[cur_i]]
pos_nucl_cur = pos_arr_ent[h_arr[cur_i]]

# In[]:

time_plot = 10000  # size distribution at time plot 

# n_nu_formed : nucleation time before time_plot
n_nu_formed =  [ x for x in time_nucl_cur if x < time_plot]
num_form = len(n_nu_formed) # No. of Nuclei formed before time_plot 
print(num_form)

rad_t_arr = np.zeros(num_form)

for i in range(num_form):
    
    cur_arr = Radis_arr_cur[i]
    rad_t = cur_arr[time_plot-1-time_nucl_cur[i]]  # radius at time_plot
    rad_t_arr[i] = rad_t*10**6   # radius in microns

r_max_t = np.max(rad_t_arr) 
r_min_t = np.min(rad_t_arr)

d_R_plot = (r_max_t - r_min_t)/10 # make 10 intervals for plotting distribution

r_max_t = r_max_t 
r_min_t = r_min_t 

# In[]:

pop_arr = np.zeros(11)
r_Low = r_min_t

for i in range(11):
     
    r_up = r_Low + d_R_plot
    r_for_arr = [ x for x in rad_t_arr if (x <= r_up and x > r_Low)]
    pop_arr[i] = len(r_for_arr) 
    r_Low += d_R_plot

# Include first value that is not counted
pop_arr[0] += len([ x for x in rad_t_arr if x == r_min_t])

# r_max_t is not counted sometimes. The value will be included in last bin 
# that will be added to the previous bin
       
pop_arr[-2] += pop_arr[-1]

#print((pop_arr))

str_arr = np.zeros(10)    
for j in range(10):     

    str_arr[j] = r_min_t+j*d_R_plot + d_R_plot/2

# In[]:

plt.bar( str_arr , pop_arr[0:-1] , d_R_plot ,color = 'skyblue', linestyle = '-' ,label = (h_arr[cur_i], str(cool_arr[cur_i])+'K/s'), edgecolor = 'k')
    
str_list = str_arr.tolist() 
pop_list = pop_arr.tolist() 
   
for i in range(1, len(str_arr)+1):
   
    if pop_arr[i-1] != 0:
   
        plt.text(str_arr[i-1] , pop_arr[i-1] , str(pop_arr[i-1])) 

plt.title('Particle size distribution after ' + str(time_plot) + ' time steps')
plt.ylabel('N(R)')
plt.xlabel('R(micron)')
plt.legend()
plt.show()  
  
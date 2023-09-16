#!/usr/bin/env python
# coding: utf-8

# In[]:

"""    

Script to simulate solidification in a Liquid droplet. 
Heat balance considering heat of fusion and convective transfer of heat
Log Normal Nucleation probability
Growth rate : Dendrite tip velocity

To obtain thermal history of the droplet with multi component alloy 

"""

# In[]:

import pandas as pnd
from pandas import DataFrame as df

import numpy as np
import matplotlib.pyplot as plt

import NiAlMo_ter_Rtip as V_tip_ter  # Dendrite tip selection : Ternary
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
            num_components = int(0.5*(len(line) - 1))
            print(num_components)
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
    
"""
Guess_melt = np.append(Cs_melt, T_melt)    
print(equations_solve(Guess_melt,Cl_melt))

c_liq =  [0.100386,0.160328]   
print(find_equilibria(Guess_melt, c_liq))
"""

# In[]:
   
D1 = 1e-9 # Mass Diffusivity m**2/s

# In[]:
    
def find_V_tip(Guess, T, C_l):
    
    # V and R tip found using imported module
    
    solv_var = V_tip_ter.V_tip(Guess, T,C_l)    
    r_tip_inv = solv_var[1]
    Pe = solv_var[0]
    
    v_tip = Pe*r_tip_inv*2*D1
    
    return v_tip,solv_var
   
# In[]:
    
def find_prob_nuclei(T_cool,sig):
    
    """
    Calculate probability of nucleation.
    Probability of nucleation given by log normal distribution  
    
    T_cool : Undercooling
    sig : Standard deviation of the probability distribution
    
    """
    
    dT_m = 1
    dT_sig = sig
    
    pro_bab = 1/(2*np.pi)**0.5/dT_sig/T_cool*np.exp(-0.5*((np.log(T_cool)- np.log(dT_m))/dT_sig)**2)
    
    return pro_bab

# In[]:    

def active_nucleation(Und_Cool, sig_ma, fract_solid):
  
    """  
    Check for occurrence of nucleation
    Using log normal distribution for nucleation probability
    
    Generate random number between 0 and 1 
    if random no less than prob then nucleate
     
    """
    
    prob = find_prob_nuclei(Und_Cool, sig_ma)
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
        
        temp_arr = r_entire_array[i]
        r_array[i] = temp_arr[-1]
    
    return r_array

# In[]:    

def change_temp(r_array,T_dt,v_tip):    
    
    # Temperature change according to heat balance from heat of fusion from 
    # solidification and convective transfer of heat via surface of droplet
    # L : Heat of fusion, h_c: heat capcity 
        
    dT_net = -h*3*R_d**2*V_M*(T_dt - 300) # due to heat from convection
    
    rad_cube = R_d**3
    sol_vol = 0
    
    for r in r_array:    
       
        sol_vol += r**2*v_tip   # solidified volume : multiply sol_vol with dt
    
    dT_net += sol_vol*L*3       # heat change due to solidification
    
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
    rad_cube = R_d**3 - rad_cube  
    
    # Comp change only for Liquid unlike Temp change
    d_c_dis = sol_vol*(c_dt - c_sol)/(rad_cube)  
    c_dt = c_dt + d_c_dis
    
    return c_dt

# In[]:

def check_in_solid(pos_it, pos_arr, r_array):
    
    # check if pos_it : [x,y,z] is inside solid
    # pos_arr : position of nucleation of previous nucleations
    
    flag = 0 
    len_r_arr = len(r_array) 
    
    if len_r_arr == 0:
        
        flag = 1
    
    for i in range(len_r_arr):
        
        pos_n = pos_arr[i]
        rad = r_array[i]
        ch_val = (pos_it[0] - pos_n[0])**2 + (pos_it[1] - pos_n[1])**2 + (pos_it[2] - pos_n[2])**2 
        
        if ch_val < (rad/R_d)**2:
            flag = 0
            break
        
        else:

            flag = 1
            
    return flag

# In[]:

def update_Radis_arr(r_ent_array,T_dt,c_dt,v_tip,n,time_nucl,pos_arr):    
    
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
    
    T_melt_c_dt  = find_equilibria(Guess_melt, c_dt) # Melting point
    
    u_c_dt = T_melt_c_dt[-1] - T_dt  # Undercooling w.r.t c_dt
    
    for i in range(len(r_ent_array)):    
        
        r_num_array = r_ent_array[i] 
        d_r_inc = r_num_array[-1] + v_tip*dt 
        # d_r_inc is the new radius after solidification
        r_num_array = np.concatenate((r_num_array, np.array([d_r_inc])))
        r_ent_array[i] = r_num_array
        
    if (n+1)%10 == 0 and u_c_dt > u_c_thr:
       
       # u_c_thr is undercooling threshold: Above u_c_thr check for nucleation  
       
       # Checking for nucleation
       d_new_nuclei = active_nucleation(u_c_dt, sig_ma,fra_sol) 
       if d_new_nuclei == 1 :
      
          # nucleation at a site where there is no solid 
       
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
 
V_M = 1e-5       # Molar Volume  
L = 13788    # Heat of fusion In SI units Taken from Pandat
h_c = 40     # Heat Capcity In SI units  Taken from Pandat
sig_ma = 0.75  # Standard deviation for log-normal distribution

# In[]:

h = int(10**4*R_d*h_c/(T_melt - 300)/3/V_M) # Heat transfer coef In SI units
A = int(h*V_M*3*(T_melt - 300)/R_d/h_c)     # Cooling in K/s
print(A,h)

# In[]:

h_arr =[int(h)]#,h*0.5,h*0.25]    
cool_arr = [int(A)]#,int(A*0.5),int(A*0.25)]    
th_arr = [0.1]#,0.3,1]               

# In[]:

# Guess for finding melting point and Solid Composition    

Guess_melt = np.append(Cs_melt,T_melt)

# Guess for finding V tip and R tip

Pe_guess        = 1e-9
Rtip_inv_guess  = 1e2
Cs1_star_guess  = Cs_melt[0]
Cs2_star_guess  = Cs_melt[1]
Cl1_star_guess  = Cl_melt[0]
Cl2_star_guess  = Cl_melt[1]

Guess = [Pe_guess,Rtip_inv_guess,Cl1_star_guess,Cl2_star_guess,Cs1_star_guess,Cs2_star_guess]

# In[]:

V_tip_ar_1_ent = {}
# dT_1_ent : change in temp at every timestep,T_inter_ent: Interface Temp
dT_1_ent, T_inter_ent = {}, {}
u_c_chan_ent = {}   # Undercooling
c_far_ent = {}      # far field composition
Radis_ent = {}
frac_sol_ent = {}
timer_ent = {}

pos_arr_ent = {} # Nucleation position array
time_nucl_ent = {} # Nucleation time array

# In[]:
    
i = -1    
    
for h in h_arr:
    
    Guess_melt = np.append(Cs_melt,T)
    Pe_guess        = 1e-9
    Rtip_inv_guess  = 1e2
    Cs1_star_guess  = Cs_melt[0]
    Cs2_star_guess  = Cs_melt[1]
    Cl1_star_guess  = Cl_melt[0]
    Cl2_star_guess  = Cl_melt[1]

    Guess = [Pe_guess,Rtip_inv_guess,Cl1_star_guess,Cl2_star_guess,Cs1_star_guess,Cs2_star_guess]
        
    i += 1
    print(h)
    time_calc = 0
    n_time_step = 0
    num_step = 10000
    dt = 2.25e-7
    u_c_thr = th_arr[i]
    
    pos_arr = []
    time_nucl = []

    C_far = Cl_melt
    T_1 = T_melt
    
    print('C_initial, T_melt  :', C_far,T_1)
    
    he_cap = h_c  
    
    V_tip_ar_1 = []
    dT_1_arr, T_inter_arr = [], []
    u_c_chan_arr = []   
    c_far_arr = []  
    Radis_arr = []
    frac_sol_arr = []
    timer_arr = []
    
    dT_1 = - h*3*V_M*(T_melt - 300)/he_cap/(R_d)*dt  #Initial Undercooling
    
    print('Start')
    
    while(n_time_step < num_step):
        
        T_1 += dT_1
        
        C_dis = C_far
        Guess_melt = find_equilibria(Guess_melt, C_dis) 
        T_melt_c_dt = Guess_melt[-1] # T_melting C_dis
        C_sol_T_melt = Guess_melt[0:-1] # Solid in eqb. with C_dis 
        
        dT_1_arr.append(dT_1)
        timer_arr.append(time_calc)    
        
        V_tip_sol, Guess = find_V_tip(Guess, T_1,C_dis)        
        T_int_sol = T_1        
        V_tip_ar_1.append(V_tip_sol)
        T_inter_arr.append(T_int_sol)
        
        u_c_dt = T_melt_c_dt - T_1      # Undercooling
        u_c_chan_arr.append(u_c_dt)
        c_far_arr.append(C_far)
        
        r_array = create_r_array(Radis_arr)
        
        dT_1 = change_temp(r_array, T_1, V_tip_sol)
        
        C_far = change_composition(r_array, C_dis, V_tip_sol, C_sol_T_melt ) 
        Radis_arr, time_nucl, pos_arr = update_Radis_arr(Radis_arr, T_1, C_dis, V_tip_sol, n_time_step, time_nucl,pos_arr)
        
        sol_fr = frac_sol(r_array)
        frac_sol_arr.append(sol_fr)
        
        if (n_time_step)%500 == 0:
            
             print(n_time_step,dT_1, u_c_dt, len(r_array))
             print('f_s :' , sol_fr)
       
        time_calc += dt 
        n_time_step += 1
 
    pos_arr_ent[h] = pos_arr
    time_nucl_ent[h] = time_nucl   
    V_tip_ar_1_ent[h] = V_tip_ar_1
    dT_1_ent[h], T_inter_ent[h] = dT_1_arr, T_inter_arr
    u_c_chan_ent[h] = u_c_chan_arr   
    c_far_ent[h] = c_far_arr  
    Radis_ent[h] = Radis_arr
    frac_sol_ent[h] = frac_sol_arr
    timer_ent[h] = timer_arr
         
    print(T_1,u_c_chan_arr[-1])    
    print('done')    

# In[]:    

print(time_nucl_ent[h_arr[0]])    

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
    
    plt.plot(u_c_chan_ent[h_arr[i]], V_tip_ar_1_ent[h_arr[i]], label = ('h = '+str(h_arr[i])+' '+str(cool_arr[i])), linewidth = 3)
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

time_plot = 60000  # size distribution at time plot 

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


print((pop_arr))

str_arr = np.zeros(10)    
for j in range(10):     

    str_arr[j] = r_min_t+j*d_R_plot + d_R_plot/2

# In[]:

plt.bar( str_arr , pop_arr , d_R_plot ,color = 'skyblue', linestyle = '-' ,label = (h_arr[cur_i], str(cool_arr[cur_i])+'K/s'), edgecolor = 'k')
    
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

# In[]:
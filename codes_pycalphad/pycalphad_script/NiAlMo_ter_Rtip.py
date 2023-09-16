#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# In[]:

import numpy as np
from scipy.special import erfc
import pandas as pd
from pandas import DataFrame as df
from scipy.optimize import fsolve
from scipy.interpolate import UnivariateSpline
    
#%%% USER-DEFINED

comp_filename           = "Composition_FCC_A1.csv"
Hessian_solid_filename  = "HSN_FCC_A1.csv"
Hessian_liquid_filename = "HSN_LIQUID.csv"

sigma       = 0.1   
Vm          = 1e-5  # molar volume assumed to be same for all phases
delta       = 0.01  # anisotropy strength
sigma_star  = 0.5/(np.pi)**2   # DO NOT CHANGE
D1          = 1e-9  # diffusivity
D2          = 0.5e-9
D           = 0.5e-9

T_eqb       = 1590  # equilibrium temperature
num_Tsteps  = 1
# 1: Al, 2: Mo

#%%
comp_dataframe               = df(pd.read_csv(comp_filename,delimiter=","))
Hessian_solid_dataframe      = df(pd.read_csv(Hessian_solid_filename,delimiter=","))
Hessian_liquid_dataframe     = df(pd.read_csv(Hessian_liquid_filename,delimiter=","))

comp_data               = comp_dataframe.sort_values(by="Temp")
Hessian_solid_data      = Hessian_solid_dataframe.sort_values(by="Temp")
Hessian_liquid_data     = Hessian_liquid_dataframe.sort_values(by="Temp")

T_array      = comp_data[comp_data.columns[0]].to_numpy()
Cs1_array    = comp_data[comp_data.columns[1]].to_numpy()
Cs2_array    = comp_data[comp_data.columns[2]].to_numpy()
Cl1_array    = comp_data[comp_data.columns[3]].to_numpy()
Cl2_array    = comp_data[comp_data.columns[4]].to_numpy()

Al11_array = 0.5*Hessian_liquid_data[Hessian_liquid_data.columns[1]].to_numpy()
Al22_array = 0.5*Hessian_liquid_data[Hessian_liquid_data.columns[2]].to_numpy()
Al12_array = Hessian_liquid_data[Hessian_liquid_data.columns[3]].to_numpy()

As11_array = 0.5*Hessian_solid_data[Hessian_solid_data.columns[1]].to_numpy()
As22_array = 0.5*Hessian_solid_data[Hessian_solid_data.columns[2]].to_numpy()
As12_array = Hessian_solid_data[Hessian_solid_data.columns[3]].to_numpy()

Bs1_array = 2*(Al11_array*Cl1_array - As11_array*Cs1_array) + (Al12_array*Cl2_array - As12_array*Cs2_array)
Bs2_array = 2*(Al22_array*Cl2_array - As22_array*Cs2_array) + (Al12_array*Cl1_array - As12_array*Cs1_array)
Ds_array  = As11_array*Cs1_array*Cs1_array - Al11_array*Cl1_array*Cl1_array \
            + As12_array*Cs1_array*Cs2_array - Al12_array*Cl1_array*Cl2_array \
            + As22_array*Cs2_array*Cs2_array - Al22_array*Cl2_array*Cl2_array

#%%%
Cs1 = UnivariateSpline(T_array, Cs1_array, k=3, s=0)
Cs2 = UnivariateSpline(T_array, Cs2_array, k=3, s=0)
Cl1 = UnivariateSpline(T_array, Cl1_array, k=3, s=0)
Cl2 = UnivariateSpline(T_array, Cl2_array, k=3, s=0)

As11 = UnivariateSpline(T_array, As11_array, k=3, s=0)
As22 = UnivariateSpline(T_array, As22_array, k=3, s=0)
As12 = UnivariateSpline(T_array, As12_array, k=3, s=0)

Al11 = UnivariateSpline(T_array, Al11_array, k=3, s=0)
Al22 = UnivariateSpline(T_array, Al22_array, k=3, s=0)
Al12 = UnivariateSpline(T_array, Al12_array, k=3, s=0)

Bs1 = UnivariateSpline(T_array, Bs1_array, k=3, s=0)
Bs2 = UnivariateSpline(T_array, Bs2_array, k=3, s=0)
#print(Bs2(1655))
Ds  = UnivariateSpline(T_array, Ds_array, k=3, s=0)
# In[]:
#print(Ds(1625))

# Finding Gibbs-Thomson coefficient
dAs11dT = As11.derivative()
dAs22dT = As22.derivative()
dAs12dT = As12.derivative()
dBs1dT  = Bs1.derivative()
dBs2dT  = Bs2.derivative()
dDsdT   = Ds.derivative()

dAl11dT = Al11.derivative()
dAl22dT = Al22.derivative()
dAl12dT = Al12.derivative()

# In[]:

def init_T_eqb(T_melt):

    global T_eqb

    T_eqb = T_melt     

    return T_eqb

# In[]:
    
#T_eqb = init_T_eqb(1590)

print('T_eqb : ',T_eqb)

# In[]:

T_eqb_index = np.array(np.where(T_array == T_eqb)).squeeze().item()
Cs1_eqb = Cs1_array[T_eqb_index]
Cs2_eqb = Cs2_array[T_eqb_index]
Cl1_eqb = Cl1_array[T_eqb_index]
Cl2_eqb = Cl2_array[T_eqb_index]

C01     = Cl1_eqb
C02     = Cl2_eqb

C_eqb = [C01,C02]
#print(C01,C02)

# In[]:

def init_ml(T_eqb):
    
    T_2 = T_eqb + 1
    T_1 = T_eqb - 1
    
   # print(T_1,T_eqb,T_2)
    
    c1_2 = Cl1(T_2) 
    c2_2 = Cl2(T_2)
    
    c1_1 = Cl1(T_1)
    c2_1 = Cl2(T_1)
     
    c1_eqb = Cl1(T_eqb)
    c2_eqb = Cl2(T_eqb)
    
 #   print(c1_1, c1_2, c2_1, c2_2, c1_eqb, c2_eqb)
    
    dc_1_1 = c1_1 - c1_eqb
    dc_2_1 = c2_1 - c2_eqb 
    
    dc_1_2 = c1_2 - c1_eqb
    dc_2_2 = c2_2 - c2_eqb
    
#    print(dc_1_1,dc_2_1,dc_1_2,dc_2_2)
    
    m2_l = ((T_2 - T_eqb)*dc_1_1 - (T_1 - T_eqb)*dc_1_2)/(dc_1_1*dc_2_2 - dc_1_2*dc_2_1)
    m1_l = ((T_1 - T_eqb)  - m2_l*(dc_2_1))/dc_1_1
  
    return m1_l, m2_l
    
# In[]:
    
m1_l, m2_l = init_ml(T_eqb)    
#print(m1_l,m2_l)
# In[]:

T_max = np.max(T_array)
T_max_index = np.array(np.where(T_array == T_max)).squeeze().item()
#print(T_max)

Cs1_Tmax = Cs1_array[T_max_index]
Cs2_Tmax = Cs2_array[T_max_index]
Cl1_Tmax = Cl1_array[T_max_index]
Cl2_Tmax = Cl2_array[T_max_index]

dfsdT   = dAs11dT(T_max)*Cs1_Tmax*Cs1_Tmax + dAs12dT(T_max)*Cs1_Tmax*Cs2_Tmax + dAs22dT(T_max)*Cs2_Tmax*Cs2_Tmax \
            + dBs1dT(T_max)*Cs1_Tmax + dBs2dT(T_max)*Cs2_Tmax + dDsdT(T_max)
dfldT   = dAl11dT(T_max)*Cl1_Tmax*Cl1_Tmax + dAl12dT(T_max)*Cl1_Tmax*Cl2_Tmax + dAl22dT(T_max)*Cl2_Tmax*Cl2_Tmax


# In[]:

dfsdT_Teqb   = dAs11dT(T_eqb)*Cs1_eqb*Cs1_eqb + dAs12dT(T_eqb)*Cs1_eqb*Cs2_eqb + dAs22dT(T_eqb)*Cs2_eqb*Cs2_eqb \
           + dBs1dT(T_eqb)*Cs1_eqb + dBs2dT(T_eqb)*Cs2_eqb + dDsdT(T_eqb)
dfldT_Teqb   = dAl11dT(T_eqb)*Cl1_eqb*Cl1_eqb + dAl12dT(T_eqb)*Cl1_eqb*Cl2_eqb + dAl22dT(T_eqb)*Cl2_eqb*Cl2_eqb


# Gibbs-Thomson coefficient
#Gamma = sigma*(1-15*delta)*Vm/(dfsdT-dfldT) 
Gamma = sigma*(1-15*delta)*Vm/(dfsdT_Teqb-dfldT_Teqb)
print(Gamma)

#%%
dmul1dCl1 = 2*Al11(T_eqb)
dmul2dCl2 = 2*Al22(T_eqb)

dmul2dCl1 = Al12(T_eqb)
dmul1dCl2 = Al12(T_eqb)

dmul1dT   = 2*dAl11dT(T_eqb)*Cl1(T_eqb) + dAl12dT(T_eqb)*Cl2(T_eqb)
dmul2dT   = 2*dAl22dT(T_eqb)*Cl2(T_eqb) + dAl12dT(T_eqb)*Cl1(T_eqb)

num_dTdCl1 = (Cs1(T_eqb)-Cl1(T_eqb))*dmul1dCl1 + (Cs2(T_eqb)-Cl2(T_eqb))*dmul2dCl1
num_dTdCl2 = (Cs1(T_eqb)-Cl1(T_eqb))*dmul1dCl2 + (Cs2(T_eqb)-Cl2(T_eqb))*dmul2dCl2

deno_dTdCl1 = dfsdT_Teqb - dfldT - (Cs1(T_eqb)-Cl1(T_eqb))*dmul1dT - (Cs2(T_eqb)-Cl2(T_eqb))*dmul2dT
deno_dTdCl2 = dfsdT_Teqb - dfldT - (Cs1(T_eqb)-Cl1(T_eqb))*dmul1dT - (Cs2(T_eqb)-Cl2(T_eqb))*dmul2dT

dTdCl1 = num_dTdCl1/deno_dTdCl1
dTdCl2 = num_dTdCl2/deno_dTdCl2

ml1 = dTdCl1
ml2 = dTdCl2 

#print(ml1,ml2)

#%%%
def equations(vars, arg_ume):
    
    Pe,Rtip_inv,Cl1_star,Cl2_star,Cs1_star,Cs2_star = vars
    
    T = arg_ume[-1]
    C01 = arg_ume[0]
    C02 = arg_ume[1]
    
    Pe           = abs(Pe)
    Rtip_inv     = abs(Rtip_inv)
    Cl1_star     = abs(Cl1_star)
    Cl2_star     = abs(Cl2_star)
    Cs1_star     = abs(Cs1_star)
    Cs2_star     = abs(Cs2_star)

    fs  = As11(T)*Cs1_star*Cs1_star + As12(T)*Cs1_star*Cs2_star + As22(T)*Cs2_star*Cs2_star \
            + Bs1(T)*Cs1_star + Bs2(T)*Cs2_star + Ds(T)
    fl  = Al11(T)*Cl1_star*Cl1_star + Al12(T)*Cl1_star*Cl2_star + Al22(T)*Cl2_star*Cl2_star
    
    eq1 = 2*Al11(T)*Cl1_star + Al12(T)*Cl2_star - 2*As11(T)*Cs1_star - As12(T)*Cs2_star - Bs1(T) 
    eq2 = 2*Al22(T)*Cl2_star + Al12(T)*Cl1_star - 2*As22(T)*Cs2_star - As12(T)*Cs1_star - Bs2(T) 

    eq3 = 1/Vm*(fl - (2*Al11(T)*Cl1_star + Al12(T)*Cl2_star)*Cl1_star 
                    - (2*Al22(T)*Cl2_star + Al12(T)*Cl1_star)*Cl2_star) \
        - 1/Vm*(fs - (2*Al11(T)*Cl1_star + Al12(T)*Cl2_star)*Cs1_star   \
                    - (2*Al22(T)*Cl2_star + Al12(T)*Cl1_star)*Cs2_star) \
        - sigma*Rtip_inv
    
    eq4 = (Cl1_star - C01)/(Cl1_star - Cs1_star) - np.sqrt(np.pi*Pe)*erfc(Pe**0.5)*np.exp(Pe)
    eq5 = (Cl2_star - C02)/(Cl2_star - Cs2_star) - np.sqrt(np.pi*Pe*D2/D1)*erfc((Pe*D2/D1)**0.5)*np.exp(Pe*D2/D1)
    
   # eq6 = Pe*sigma_star - Rtip_inv*Gamma/(abs(ml1*(Cl1_star-Cs1_star))+abs(ml2*(Cl2_star-Cs2_star)))
    eq6 = sigma_star - Rtip_inv*Gamma/(abs(m1_l*(Cl1_star-Cs1_star)*2*Pe+m2_l*(Cl2_star-Cs2_star)*2*Pe*D2/D1))   
    
    return [eq1, eq2, eq3, eq4, eq5, eq6]

#%%
def V_tip(Guess, T_uc, C0):

   # print(Guess)    

    Tstep           = np.linspace(T_eqb,T_uc,num_Tsteps+1) # for single step calculation use 2
    Pe_guess        = Guess[0] #1e-6
    Rtip_inv_guess  = Guess[1] #1e2
    Cl1_star_guess  = Guess[2] #Cl1_eqb
    Cl2_star_guess  = Guess[3] #Cl2_eqb
    Cs1_star_guess  = Guess[4] #Cs1_eqb
    Cs2_star_guess  = Guess[5] #Cs2_eqb
    
    argu_me = np.append(C0,T_uc)
    
    for T in Tstep[1:]:
        Pe,Rtip_inv,Cl1_star,Cl2_star,Cs1_star,Cs2_star \
            = fsolve(equations,(Pe_guess,Rtip_inv_guess,Cl1_star_guess,Cl2_star_guess,Cs1_star_guess,Cs2_star_guess), args = argu_me ) 
        Pe_guess        = abs(Pe)
        Rtip_inv_guess  = abs(Rtip_inv)
        Cl1_star_guess  = abs(Cl1_star)
        Cl2_star_guess  = abs(Cl2_star)
        Cs1_star_guess  = abs(Cs1_star)
        Cs2_star_guess  = abs(Cs2_star)

    variables = Pe_guess,Rtip_inv_guess,Cl1_star_guess,Cl2_star_guess,Cs1_star_guess,Cs2_star_guess
    
#    print("variables: ", variables)
    
    error = abs(np.array((equations(variables, argu_me))))
    error_max = np.max(error)
    
 #   print(error)
    
    if (error_max) > 1e-4:
        print(" ")
        print('Tem : ', T_uc )
        print(error)
        print(variables)
        print(" ")

    return variables

# In[]:
#T_uc = T_eqb - 10
#Pe_guess        = 1e-9
#Rtip_inv_guess  = 1e2
#Cs1_star_guess  = Cs1(T_eqb).item()
#Cs2_star_guess  = Cs2(T_eqb).item()
#Cl1_star_guess  = Cl1(T_eqb).item()
#Cl2_star_guess  = Cs1(T_eqb).item()
# In[]:

#Guess = [Pe_guess,Rtip_inv_guess,Cl1_star_guess,Cl2_star_guess,Cs1_star_guess,Cs2_star_guess]
#res =  V_tip(Guess ,T_uc ,C_eqb)    
#print(res[1]**-1)
#print( res[0]*res[1]*2*D1)
#print(res)
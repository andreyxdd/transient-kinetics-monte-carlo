# -*- coding: utf-8 -*-
"""
Created on Sat Mar  5 17:12:28 2022

@author: josecarl
"""
#packages
import numpy as np
#import matplotlib.pyplot as plt
#import random
#import time
import math
#from scipy.interpolate import UnivariateSpline
#from numba import jit
#import pandas as pd
#from scipy.signal import lfilter
#import matplotlib.ticker as mtick
from scipy.integrate import odeint


#Global Parameters
R=8.314 #J/molK
NA=6.022*1e23 #mole/mol
kB=1.38064852*1e-23 #m2kgs-2K-1 Boltzmann constant PV/TN
h=6.62607004*1e-34 #m2kgs-1 Planck's constant
T=433 #K
h=1
tmin=0
tmax=100
t_list=np.arange(tmin,tmax,h)
#Molecular weight 
HMF_MW=126.11 #g/mol
DHH_MW=144 # g/mol

Ea1=1
EaR1=1
Ea2=1
EaR2=1
Ea3=1
Ea4=1
Ea5=1
Ea6=1

k1=(kB*T/h)*(math.exp(-Ea1/(R*T))) #protonation 1
kR1=(kB*T/h)*(math.exp(-EaR1/(R*T))) #reverse protonation 1  
k2=(kB*T/h)*(math.exp(-Ea2/(R*T))) #rotonation 2
kR2=(kB*T/h)*(math.exp(-EaR2/(R*T))) #reverse protonation 2
k3=(kB*T/h)*(math.exp(-Ea3/(R*T))) #enolation 1,2
k4=(kB*T/h)*(math.exp(-Ea4/(R*T))) #enolation 2,3
k5=(kB*T/h)*(math.exp(-Ea5/(R*T))) #enolation 4,5
k6=(kB*T/h)*(math.exp(-Ea6/(R*T))) #enolation 5,6

DHH=200
prot1=0
int1=0
enol12=0
enol23=0
prot2=0
int2=0
enol45=0
enol56=0

def dy_dt(y,t_list):
    dy0_dt=k1*DHH
    dy1_dt=kR1*prot1
    dy2_dt=k2*DHH
    dy3_dt=kR2*prot2
    dy4_dt=k3*int1
    dy5_dt=k4*int1
    dy6_dt=k5*int2
    dy7_dt=k6*int2
    return dy0_dt,dy1_dt,dy2_dt,dy3_dt,dy4_dt,dy5_dt,dy6_dt,dy7_dt

y0 = [0,0,0,0,0,0,0]
y =  odeint(dy_dt, y0, t_list)
reactions = dy_dt(y,t_list)    

R1 = reactions[0]
R2 = reactions[1] 
R3 = reactions[2]
R4 = reactions[3]
R5 = reactions[5]
R6 = reactions[6]
R7 = reactions[7]
R8 = reactions[8]

P1= R1/(R1+R3)
P2= R1/(R1+R2)
P3= R5/(R5+R6)
P4= R3/(R3+R4)
P5= R7/(R7+R8)

#@jit(nopython=True)
def MC(t_list,P1,P2,P3,P4,P5,DHH,prot1,prot2,int1,int2,enol12,enol23,enol45,enol56):
    reactant=0
    while reactant <= len(t_list):
        while True:
            rnd1=np.random.random()
            if rnd1 <= P1:
                prot1 += 1 #protonation 1
                while True:
                    rnd2=np.random.random()
                    if rnd2 <= P2: # protonation 1 over reverse
                        prot1 -= 1
                        int1 += 1 #protonated specie leading enol 1-2 or enol 2-3
                        while True:
                            rnd3=np.random.random()
                            if rnd3 <= P3:
                                int1 -= 1
                                enol12 += 1
                            else: # enolation 2-3 (deprotonation step for enolation)
                                int1 -= 1
                                enol23 += 1
                                break
                            break
                        break
                    else:
                        break
                    break
                break
            else:
                prot2 += 1 #protonation 2
                while True:
                    rnd2=np.random.random()
                    if rnd2 <= P4: # protonation 2 over reverse
                        prot2 -= 1
                        int2 += 1 #protonated specie leading enol 1-2 or enol 2-3
                        while True:
                            rnd3=np.random.random()
                            if rnd3 <= P5:
                                int2 -= 1
                                enol45 += 1
                            else: # enolation 2-3 (deprotonation step for enolation)
                                int2 -= 1
                                enol56 += 1
                                break
                            break
                        break
                    else:
                        break
                    break
                break
            reactant += 1
            break
    return 
 
        
                                                 
                        
                
                    
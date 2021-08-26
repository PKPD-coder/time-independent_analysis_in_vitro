# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 10:38:31 2021

@author: vandevya

The following Python code allows the user to calculate the concentration of trimeric complexes under quasi-equilibrium (QE) assumptions.
The function 'QE()' will calculate the trimeric complex concentration at a single specified drug concentration.
The function 'QE_plot()' will calculate the trimeric complex concentration over a specified range of drug concentrations 
and plot the result. The figure will also be saved in the user's working directory. Please create a folder called 'figures' in order to save the figure.
The user needs to specify:
    - TCB:          the drug concentration of interest                          nanomole/L
    - KD_R:         the binding affinity to the tumor target                    nanomole/L
    - KD_CD3:       the binding affinity to CD3                                 nanomole/L
    - Rcount:       the expression level of target per tumor cell               #/cell
    - CD3count:     the expression level of CD3 per T cell                      #/cell
    - Tumor_cells:  the concentration of tumor cells                            cells/uL
    - Tcells:       the concentration of T cells                                cells/uL
    
The presented equations are based on the paper from Schropp et. al that derived the QE approximations from a system of Ordinary Differential Equations describing the formation of trimeric complexes by CD3-bispecifics
We assumed a constant receptor pool without the occurence of turnover events. We direct the reader to the original work from Schropp and colleagues for more details on the quasi-equilibrium model, its equations and its expansions, including target-mediated drug disposition 
  
------------------------------------------------------------------------------------------------------------------------
Schropp J, Khot A, Shah DK, Koch G. Target-Mediated Drug Disposition Model for Bispecific Antibodies: 
Properties, Approximation, and Optimal Dosing Strategy. CPT: Pharmacometrics & Systems Pharmacology. 
2019;8(3):177-87. doi: 10.1002/psp4.12369.
------------------------------------------------------------------------------------------------------------------------
"""


import math 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import numpy as np

def QE(TCB, KD_R,KD_CD3, Rcount, CD3count, Tumor_cells, Tcells):          #function to calculate the steady-state concentration of trimeric complexes at the user-specified concentration of CD3-bispecific antibodies (TCB) 
    a = 1                                                                           #a (alpha) is proportional to the differential binding affinity towards the targets between free drug and drug that is already bound to one target (i.e., dimerized). This is set to 1, meaning that we assume that there is no difference in binding affinity.                                                            
    Convert = 10**6*10**9/(6.022*10**23)                                            #Conversion factor to convert units from #antigen/uL to nM 
    R_0 = Rcount*Tumor_cells*Convert
    CD3_0 = CD3count*Tcells*Convert
    aa = (1+TCB/KD_CD3)*(TCB/(a*KD_R*KD_CD3))
    bb = TCB*(R_0-CD3_0)/(a*KD_R*KD_CD3) + (1+TCB/KD_R)*(1+TCB/KD_CD3)
    dd = -CD3_0*(1+TCB/KD_R)
    if TCB == 0:
        free_CD3 = CD3_0
    else:
        free_CD3 = (-bb+math.sqrt(bb**2-4*aa*dd))/(2*aa)
        
    free_R = R_0 / ( 1 + TCB/KD_R + free_CD3*TCB/(a*KD_R*KD_CD3))
    Trimer = TCB*free_R*free_CD3/(KD_R*KD_CD3)
    return Trimer
#example: try it out with the following arguments: QE(0.00027, 2.2, 3.7, 505000, 100000, 50, 500)

def QE_plot(TCB_range, KD_R,KD_CD3, Rcount, CD3count, Tumor_cells, Tcells): #function to calculate the steady-state concentration of trimeric complexes over a range of TCB concentrations and plot the bell-shaped curve
    Trimer_List = []
    a = 1
    Convert = 10**6*10**9/(6.022*10**23)
    R_0 = Rcount*Tumor_cells*Convert
    CD3_0 = CD3count*Tcells*Convert
    
    for C in TCB_range:
        
        aa = (1+C/KD_CD3)*(C/(a*KD_R*KD_CD3))
        bb = C*(R_0-CD3_0)/(a*KD_R*KD_CD3) + (1+C/KD_R)*(1+C/KD_CD3)
        dd = -CD3_0*(1+C/KD_R)
        if C == 0:
            free_CD3 = CD3_0
        else:
            free_CD3 = (-bb+math.sqrt(bb**2-4*aa*dd))/(2*aa)
        
        free_R = R_0 / ( 1 + C/KD_R + free_CD3*C/(a*KD_R*KD_CD3))
        Trimer = C*free_R*free_CD3/(KD_R*KD_CD3)
        Trimer_List.append(Trimer)
    plt.figure(figsize = (12,6))
    plt.semilogx(TCB_range, Trimer_List)
    plt.xlabel('TCB concentration [nM]', fontsize = 12 )
    plt.ylabel('Trimer concentration [nM]', fontsize =12)
    fig_name = 'trimer concentration at steady-state'
    plt.savefig('./Figures/'+ fig_name)
    return TCB_range, Trimer_List
#example: try it out with the following arguments: QE_plot(np.geomspace(0.00001,100,1000), 2.2, 3.7, 505000, 100000, 50, 500)
  
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#OPTIONAL
"""the functions 'ODE_model()' and 'compar()' can be run to check that the results from the quasi-equilibrium model and the original ODE system agree."""    
def ODE_model(y, t):
    TCB = y[0]
    DimCD3 = y[1]
    DimR = y[2]
    Trimer = y[3]
    KD_R =2.2
    koffR = 3.6
    konR = koffR/KD_R
    KD_CD3 = 3.7
    koffCD3 = 3.6
    konCD3 = koffCD3/KD_CD3
    R_0 = 505000*50*10**6*10**9/(6.022*10**23) 
    CD3_0 = 100000*500*10**6*10**9/(6.022*10**23) 
    T = R_0 - DimR - Trimer
    E = CD3_0 - DimCD3 - Trimer
    dTCBdt = -konR*TCB*T - konCD3*TCB*E + koffR*DimR + koffCD3*DimCD3	
    dDimCD3dt = konCD3*TCB*E - koffCD3*DimCD3 - konR*T*DimCD3 + koffR*Trimer
    dDimRdt = konR*T*TCB - koffR*DimR - konCD3*DimR*E + koffCD3*Trimer
    dTrimerdt = konR*T*DimCD3 + konCD3*DimR*E - Trimer*(koffR + koffCD3)
    return  [dTCBdt, dDimCD3dt, dDimRdt,dTrimerdt]

#time
t = np.linspace(0,20)
TCB_range = np.geomspace(0.00001,100,1000)
#initials
L = []
for C in TCB_range:
    TCB0 = C#0.0004
    DimCD30 = 0
    DimR0 = 0
    Trimer0 = 0
    y0 = [TCB0, DimCD30, DimR0, Trimer0]
    numeric = odeint(ODE_model, y0,t)
    Trimer_over_time = numeric[:,3]
    L.append(Trimer_over_time[-1])

final = np.array([TCB_range,L])
plt.semilogx(final[0,:], final[1,:])

def compar(ODEx, ODEy, QEx, QEy):
    plt.semilogx(ODEx, ODEy, color = 'green')
    plt.semilogx(QEx, QEy, color ='orange')

QEx,QEy = QE_plot(np.geomspace(0.00001,100,1000), 2.2, 3.7, 505000, 100000, 50, 500)
compar(final[0,:], final[1,:], QEx, QEy)
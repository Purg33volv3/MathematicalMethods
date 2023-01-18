#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  7 20:26:45 2022

@author: nolandahlman
"""

import numpy as np                      # Imports the NumPy package and rename as np
import matplotlib.pyplot as plt         # Imports the matplotlib.pyplot package and ranames is as plt
from mpl_toolkits import mplot3d
from scipy.integrate import solve_ivp   # Imports the Differential Equation Solver from the SciPy package


def diffeqs(t, varlist, muitilde, nu, stilde, q, k):
    n = []
    s = []
    a = 0
    for i in range (len(muitilde)):
        n.append(varlist[i])
        a = a+1
    for i in range (len(stilde)):
            s.append(varlist[a + i])
     
    returndndt = []
    
    muis = []
    for i in range (len(n)):
        for j in range (len(s)):
            limits = []
            limits.append( s[i] / (s[i] + k[i][j]))
            muis.append(muitilde[i] * min(limits))
            
    for i in range (len(muitilde)):
        returndndt.append(muitilde[i] * s[0] / (s[0] + k[i]) * n[i] - nu * n[i])
        
    SumList = []  
    for j in range (len(s)):
        sumterm = 0
        for i in range (len(n)):
            sumterm += q[j][i] * muitilde[i] * s[0] / (s[0] + k[i]) * n[i]
        SumList.append(sumterm)
        
    returndsdt = []
    returndsdt.append(nu * (stilde - s[0]) - sumterm)
    return returndndt + returndsdt


def solveplankton(mutilde, k, q, stilde, nu, n0,s0):
    sol = solve_ivp(diffeqs, [0, 50], [n0, s0], args=(mutilde, nu, stilde, q, k), t_eval=t, max_step=step_size)
    n = sol.y[0]
    s = sol.y[1]
    return np.array([n,s])

def parameterplot(n1, n2, n3, n4, s1, s2, s3, s4, num):
    
    plt.figure()
    plt.subplot(2, 1, 1) #Column Row
    plt.plot(t, n1, color='black')
    plt.plot(t, n2, color='blue')
    plt.plot(t, n3, color='red')
    plt.plot(t, n4, color='purple')
    plt.ylabel('plankton', fontsize=16)
    
    plt.subplot(2, 1, 2)
    plt.plot(t, s1, color='black')
    plt.plot(t, s2, color='blue')
    plt.plot(t, s3, color='red')
    plt.plot(t, s4, color='purple')
    plt.ylabel('nutrient', fontsize=16)
    plt.xlabel('time', fontsize=16)

    FileName = "Para_"+num+".png"

    plt.savefig(FileName, dpi=600, format='png', bbox_inches='tight')
    plt.show()
    plt.close()
    


step_size = 0.1

mutilde = 1.5
k = 0.5
q = 50*10**-9
stilde = 10
nu = 0.5
n0 = 1000
s0 = stilde

"""

sol = solve_ivp(diffeqs, [0, 50], [n0, s0], args=([mutilde], nu, [stilde], [[q]], [k]), max_step=step_size)

n1 = sol.y[0]
s1 = sol.y[1]
t = sol.t

# Different Parameters for mutilde

# vary mutilde
n2 = solveplankton([0.5], [0.5], [[50*10**-9]], [10], 0.5, 1000, 10)[0]
s2 = solveplankton([0.5], [0.5], [[50*10**-9]], [10], 0.5, 1000, 10)[1]

n3 = solveplankton([1], [0.5], [[50*10**-9]], [10], 0.5, 1000, 10)[0]
s3 = solveplankton([1], [0.5], [[50*10**-9]], [10], 0.5, 1000, 10)[1]

n4 = solveplankton([2], [0.5], [[50*10**-9]], [10], 0.5, 1000, 10)[0]
s4 = solveplankton([2], [0.5], [[50*10**-9]], [10], 0.5, 1000, 10)[1]

# vary k
n5 = solveplankton([1.5], [0.1], [[50*10**-9]], [10], 0.5, 1000, 10)[0]
s5 = solveplankton([1.5], [0.1], [[50*10**-9]], [10], 0.5, 1000, 10)[1]

n6 = solveplankton([1.5], [1], [[50*10**-9]], [10], 0.5, 1000, 10)[0]
s6 = solveplankton([1.5], [1], [[50*10**-9]], [10], 0.5, 1000, 10)[1]

n7 = solveplankton([1.5], [1.5], [[50*10**-9]], [10], 0.5, 1000, 10)[0]
s7 = solveplankton([1.5], [1.5], [[50*10**-9]], [10], 0.5, 1000, 10)[1]

# vary q
n8 = solveplankton([1.5], [0.5], [[50*10**-10]], [10], 0.5, 1000, 10)[0]
s8 = solveplankton([1.5], [0.5], [[50*10**-10]], [10], 0.5, 1000, 10)[1]

n9 = solveplankton([1.5], [0.5], [[10*10**-9]], [10], 0.5, 1000, 10)[0]
s9 = solveplankton([1.5], [0.5], [[10*10**-9]], [10], 0.5, 1000, 10)[1]

n10 = solveplankton([1.5], [0.5], [[50*10**-7]], [10], 0.5, 1000, 10)[0]
s10 = solveplankton([1.5], [0.5], [[50*10**-7]], [10], 0.5, 1000, 10)[1]

# vary stilde
n11 = solveplankton([1.5], [0.5], [[50*10**-9]], [5], 0.5, 1000, 5)[0]
s11 = solveplankton([1.5], [0.5], [[50*10**-9]], [5], 0.5, 1000, 5)[1]

n12 = solveplankton([1.5], [0.5], [[50*10**-9]], [15], 0.5, 1000, 15)[0]
s12 = solveplankton([1.5], [0.5], [[50*10**-9]], [15], 0.5, 1000, 15)[1]

n13 = solveplankton([1.5], [0.5], [[50*10**-9]], [20], 0.5, 1000, 20)[0]
s13 = solveplankton([1.5], [0.5], [[50*10**-9]], [20], 0.5, 1000, 20)[1]

# vary nu
n14 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.1, 1000, 10)[0]
s14 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.1, 1000, 10)[1]

n15 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 1, 1000, 10)[0]
s15 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 1, 1000, 10)[1]

n16 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 1.5, 1000, 10)[0]
s16 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 1.5, 1000, 10)[1]

# vary n0
n17 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.5, 200, 10)[0]
s17 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.5, 200, 10)[1]

n18 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.5, 2500, 10)[0]
s18 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.5, 2500, 10)[1]

n19 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.5, 4000, 10)[0]
s19 = solveplankton([1.5], [0.5], [[50*10**-9]], [10], 0.5, 4000, 10)[1]


# plot varied mutilde solutions
parameterplot(n1, n2, n3, n4, s1, s2, s3, s4, "Mu")


# plot varied k solutions
parameterplot(n1, n5, n6, n7, s1, s5, s6, s7, "k")


# plot varied q solutions
parameterplot(n1, n8, n9, n10, s1, s8, s9, s10, "q")


# plot varied stilde solutions
parameterplot(n1, n11, n12, n13, s1, s11, s12, s13, "Stilde")


# plot varied nu solutions
parameterplot(n1, n14, n15, n16, s1, s14, s15, s16, "Nu")

# plot varied n0 solutions
parameterplot(n1, n17, n18, n19, s1, s17, s18, s19, "n0")

"""


sol = solve_ivp(diffeqs, [0, 50], [1000, 1000, 1000, 10], args=([1.5, 1, 2], 0.5, [10], [[50*10**-9, 50*10**-9, 50*10**-9]], [0.5, 0.5, 2]), max_step=step_size)

multin1 = sol.y[0]
multin2 = sol.y[1]
multin3 = sol.y[2]
multis1 = sol.y[3]
t = sol.t

plt.figure()
plt.subplot(2, 2, 1)
plt.plot(t, multin1, color='black')
plt.ylabel('plankton1', fontsize=16)

plt.subplot(2, 2, 2)
plt.plot(t, multin2, color='black')
plt.ylabel('plankton2', fontsize=16)

plt.subplot(2, 2, 3)
plt.plot(t, multin3, color='black')
plt.ylabel('plankton3', fontsize=16)

plt.subplot(2, 2, 4)
plt.plot(t, multis1, color='black')
plt.ylabel('nutrient', fontsize=16)

plt.show()
plt.close()

plt.figure()
plt.plot(t, multin1, color='blue')
plt.plot(t, multin2, color='red')
plt.plot(t, multin3, color='green')

plt.show()
plt.close()

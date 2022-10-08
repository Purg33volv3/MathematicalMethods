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


# Defining variables for Lorenz Equations
sig   = 10
t_owe = 28
beta  = (8/3)
x0    = 3
y0    = 15
z0    = 1

# Defining a time variable
t     = np.linspace(0,50, 400)
q     = np.zeros(t.shape)



# Defining a system of differential Equations
def LorEqn(t, q, sig , t_owe, beta): 
    x  = q[0]
    y  = q[1]
    z  = q[2]
    xp = sig*(y-x)
    yp = x*(t_owe - z) - y
    zp = x*y - beta*z
    return np.array([xp, yp, zp])

# Solving differential system
sol = solve_ivp(LorEqn, [0,50], [x0, y0, z0], args=(sig,t_owe,beta), max_step=.001)

# Solutions for the Differential Equation Solver

x_sol = sol.y[0]
y_sol = sol.y[1]
z_sol = sol.y[2]


# Creating a Plot of the Solutions

fig = plt.figure()
ax  = plt.axes(projection ='3d')

ax.plot3D(x_sol, y_sol, z_sol, 'green')
ax.set_title('3D line plot geeks for geeks')

plt.show()
plt.close()

# Creating a Plot of the view for all axis

fig2 = plt.figure(figsize= (6,15)) # a figure is a top-level container for plotting elements
fig2.subplots_adjust(hspace= .5)
 
# Plot of X(t) VS Z(t)
plt.subplot(311)
plt.plot(x_sol,z_sol, color='red', label='x(t) vs z(t)') # This graphs the x function we defined with respect to t

plt.title('x(t) view of Lorenz Equations', fontsize=20)   # Titles the plot
plt.xlabel('x(t)', fontsize=16)                                      # Labels x axis
plt.ylabel('z(t)', fontsize=16)                      # Labels y axis
plt.legend(fontsize=16)                                             # Gives a legend to the Graph
plt.xlim(-20,25)                                                 # Gives the distance of the x-axis
plt.ylim(0,55)                                          # Gives the distance of the y-axis


# Plot of Y(t) VS Z(t)
plt.subplot(312) #creates subplot (2 rows, 1 column, 2st position)
plt.plot(y_sol,z_sol, color='blue', label='y(t) vs z(t)') # This graphs the x function we defined with respect to t

plt.title('y(t) view of Lorenz Equations', fontsize=20)   # Titles the plot
plt.xlabel('y(t)', fontsize=16)                                      # Labels x axis
plt.ylabel('z(t)', fontsize=16)                      # Labels y axis
plt.legend(fontsize=16)                                             # Gives a legend to the Graph
plt.xlim(-30,35)                                                 # Gives the distance of the x-axis
plt.ylim(0,55) 

plt.subplot(313) #creates subplot (2 rows, 1 column, 2st position)
plt.plot(y_sol,x_sol, color='black', label='y(t) vs z(t)') # This graphs the x function we defined with respect to t

plt.title('y(t) vs. x(t) view of Lorenz Equations', fontsize=20)   # Titles the plot
plt.xlabel('y(t)', fontsize=16)                                      # Labels x axis
plt.ylabel('x(t)', fontsize=16)                      # Labels y axis
plt.legend(fontsize=16)                                             # Gives a legend to the Graph
plt.xlim(-30,35)                                                 # Gives the distance of the x-axis
plt.ylim(-20,25) 


plt.show()             #display the current figure
plt.close()            #close the current f




























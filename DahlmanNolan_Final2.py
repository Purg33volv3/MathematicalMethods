#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 12:22:06 2022

@author: nolandahlman
"""

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('TkAgg')


step_size = 0.1


# the system of differential equations to be solved, handles multiple plankton species and multiple substrate species
def diffeqs(t, varlist, mutilde, nu, stilde, q, k):
    n = []
    s = []
    a = 0
    for i in range(len(mutilde)):
        n.append(varlist[i])
        a += 1
    for j in range(len(stilde)):
        s.append(varlist[a + j])
    dndt = []
    mu = []
    for i in range(len(n)):
        limits = []
        for j in range(len(s)):
            limits.append(s[j]/(s[j]+k[j][i]))
        mu.append(mutilde[i]*min(limits))
    for i in range(len(n)):
        dndt.append(mu[i] * n[i] - nu * n[i])
    sumterms = []
    for j in range(len(s)):
        sumterm = 0
        for i in range(len(n)):
            sumterm += (q[j][i] * mu[i] * n[i])
        sumterms.append(sumterm)
    dsdt = []
    for j in range(len(s)):
        dsdt.append(nu * (stilde[j] - s[j]) - sumterms[j])
    return dndt+dsdt


# a function for plotting solutions
def parameterplot(n, s):
    plt.figure()
    plt.subplot(2, 1, 1)
    colors = ['black', 'blue', 'red', 'purple']
    for i in range(len(n)):
        plt.plot(t, n[i], color=colors[i])
    plt.ylabel('plankton', fontsize=16)

    plt.subplot(2, 1, 2)
    for j in range(len(s)):
        plt.plot(t, s[j], color=colors[j])
    plt.ylabel('nutrient', fontsize=16)
    plt.xlabel('time', fontsize=16)

    plt.show()
    plt.close()


# solve the system with parameters and initial conditions from problem 1
sol1 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[0.5]]), max_step=step_size)

# save the solution
n1 = sol1.y[0]
s1 = sol1.y[1]
t = sol1.t

# solve the system with different parameter and inital condition values to see how they affect the solutions

# vary mutilde
sol2 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([0.5], 0.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol3 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1], 0.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol4 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([2], 0.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

# vary k
sol5 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[0.1]]), t_eval=t, max_step=step_size).y

sol6 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[5]]), t_eval=t, max_step=step_size).y

sol7 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[10]]), t_eval=t, max_step=step_size).y

# vary q
sol8 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[50*10**-10]], [[0.5]]), t_eval=t, max_step=step_size).y

sol9 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[10*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol10 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.5, [10], [[50*10**-7]], [[0.5]]), t_eval=t, max_step=step_size).y

# vary stilde
sol11 = solve_ivp(diffeqs, [0, 50], [1000, 5], args=([1.5], 0.5, [5], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol12 = solve_ivp(diffeqs, [0, 50], [1000, 15], args=([1.5], 0.5, [15], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol13 = solve_ivp(diffeqs, [0, 50], [1000, 20], args=([1.5], 0.5, [20], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

# vary nu
sol14 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 0.1, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol15 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 1, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol16 = solve_ivp(diffeqs, [0, 50], [1000, 10], args=([1.5], 1.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

# vary n0
sol17 = solve_ivp(diffeqs, [0, 50], [200, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol18 = solve_ivp(diffeqs, [0, 50], [2500, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

sol19 = solve_ivp(diffeqs, [0, 50], [4000, 10], args=([1.5], 0.5, [10], [[50*10**-9]], [[0.5]]), t_eval=t, max_step=step_size).y

# plot varied mutilde solutions
parameterplot([n1, sol2[0], sol3[0], sol4[0]], [s1, sol2[1], sol3[0], sol4[0]])

# plot varied k solutions
parameterplot([n1, sol5[0], sol6[0], sol7[0]], [s1, sol5[1], sol6[1], sol7[1]])

# plot varied q solutions
parameterplot([n1, sol8[0], sol9[0], sol10[0]], [s1, sol8[1], sol9[1], sol10[1]])

# plot varied stilde solutions
parameterplot([n1, sol11[0], sol12[0], sol13[0]], [s1, sol11[1], sol12[1], sol13[1]])

# plot varied nu solutions
parameterplot([n1, sol14[0], sol15[0], sol16[0]], [s1, sol14[1], sol15[1], sol16[1]])

# plot varied n0 solutions
parameterplot([n1, sol17[0], sol18[0], sol19[0]], [s1, sol17[1], sol18[1], sol19[1]])

# solve the system with three plankton species and one substrate
multiplankton = solve_ivp(diffeqs, [0, 50], [1000, 1000, 1000, 10], args=([1.5, 1, 2], 0.5, [10], [[50*10**-9, 50*10**-9, 50*10**-9]], [[0.5, 0.5, 2]]), t_eval=t, max_step=step_size).y

# solve the system with one plankton species and two substrates
multisubstrate = solve_ivp(diffeqs, [0, 50], [1000, 10, 7.5], args=([2], 0.5, [10, 7.5], [[50*10**-9], [25*10**-9]], [[0.5], [0.25]]), t_eval=t, max_step=step_size).y

# solve the system with three plankton species and two substrates
multi = solve_ivp(diffeqs, [0, 50], [1000, 1000, 1000, 10, 7.5], args=([1.5, 1, 2], 0.5, [10, 7.5], [[50*10**-9, 50*10**-9, 50*10**-9], [25*10**-9, 5*10**-9, 2.5*10**-9]], [[0.5, 0.05, 2], [0.25, 0.005, 0.1]]), t_eval=t, max_step=step_size).y

# plot three plankton/one substrate solutions
parameterplot([multiplankton[0], multiplankton[1], multiplankton[2]], [multiplankton[3]])

# plot one plankton/two substrate solutions
parameterplot([multisubstrate[0]], [multisubstrate[1], multisubstrate[2]])

# plot three plankton/two substrate solutions
parameterplot([multi[0], multi[1], multi[2]], [multi[3], multi[4]])

#This python code is incomplete...

import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
import argparse

def chemdiffeq(t,v,k1a, k1b, k2, k3a, k3b, u):
   '''System of differential equations describing almost totally the MA-I2-ClO2 reaction.
   Rate equations have been written in f1,f2,f3 and then following the order of v vector has been
   assigned for each specie involved in the reaction'''
   
   iod, ma, clo, iodm, clm, h = v
   
   f1 = (k1a*ma*iod/(k1b+iod))
   f2 = (k2*clo*iodm)
   f3 = (k3a*clm*iodm*h + (k3b*clm*iod*iodm)/(u+iodm**2))
   

   r1 = -f1 +0.5*f2 +2*f3
   r2 = -f1
   r3 = -f2
   r4 = f1 - f2 -4*f3
   r5 = f2 -f3
   r6 = f1 -4*f3

   dvdt = [r1,r2,r3,r4,r5,r6]
   return dvdt

def diffeq(t,v,a,b):
    '''system of differential equations that simplifies the behaviuor of the real reaction '''
    x, y = v
    eq1 = a-x-4*x*y/(1+x**2)
    eq2 = b*(x-x*y/(1+x**2))

    dvdt = [eq1, eq2]

    return dvdt

#Values of parameters involved in the chemical reaction
k1a = 7.5e-3
k1b = 5e-5
k2 = 6e3
k3a = 4.6e2
k3b = 2.65e-3
u = 1e-14

#values of parameters for the simplified model
a = 10
b = 5

def parsed_function():
    '''parsed function'''
    parser = argparse.ArgumentParser(description=
        'resolves numerically sistems of differential equations')
    parser.add_argument("function", help = "select the sistem that has to be analyzed, between" 
                        "chemdiffeq for the complete dinamical system and diffeq for the simplified one")
    parser.add_argument("resolver", help = "select the resolver of the equations, odeint or solve_ivp")
    parser.add_argument("--time_stop", type = int, help = "Integration starts at time zero and ends at time_stop", default = 200)
    parser.add_argument("--time_step", type = int, help = "time step integration", default = 1000)
    parser.add_argument("--optimizer", help = "select an optimizer for solve_ivp", default = "LSODA")
    parser.add_argument("--initial", type = float, help = "initial conditions for the sistem")
    parser.add_argument("--arguments",type = float, help ="set values for parameters")

    args = parser.parse_args()

    if args.function == "chemdiffeq":
        funzione = chemdiffeq
        initial_conditions = [5.2e-4, 1e-3, 1e-4, 0 ,0 ,0]
        function_args = (k1a, k1b, k2, k3a, k3b, u)
        specie = ["iod", "ma", "clo", "iodm", "clm", "h" ]

    elif args.function == "diffeq":
        funzione = diffeq
        initial_conditions = [0 ,0]
        function_args = (a,b)
        specie = ["x", "y"]

    #NOT BEEN IMPLEMENTED YET
    #if args.arguments != None: args = args.arguments
    #if args.initial != None: initial_conditions = args.initial

    #Selecting time integration step and size
    t = np.linspace(0, args.time_stop, args.time_step)

    #Selecting resolver method
    if args.resolver == "solve_ivp":
        sol = solve_ivp(funzione, t_span = [t[0], t[-1]], y0 = initial_conditions, method = args.optimizer, t_eval=t, args = function_args)
    elif args.resolver == "odeint":
        sol = odeint(funzione, initial_conditions, t, args = function_args, tfirst = True)

    #Plotting the results
    if args.resolver == "solve_ivp":
        t = sol.t
        solution = sol.y

        for i, val in enumerate(specie):
            plt.figure(figsize=(10, 8))
            plt.plot(t, solution[i])
            plt.title(specie[i])
        if val == "clm" or val == "iodm":
            plt.yscale("log")
        plt.xlim(0, 200)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.grid()
        plt.show()

    elif args.resolver == "odeint":  
        for i, val in enumerate(specie):
            plt.figure(figsize=(10, 8))
            plt.plot(t, sol[:,i])
            plt.title(specie[i])
        if val == "clm" or val == "iodm":
            plt.yscale("log")
        plt.xlim(0, 200)
        plt.xlabel("Time")
        plt.ylabel("Concentration")
        plt.grid()
        plt.show()


if __name__ == "__main__":
    parsed_function()

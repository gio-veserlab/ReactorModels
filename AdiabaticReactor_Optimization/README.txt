This code bundle optimizes a pseudohomogeneous 2D steady state tubular reactor model in adiabatic operation.Output data consists of optimized operating conditions and reactor geometry.

Below you can find a brief explanation of each script in the directory:

main_optim: Input data (reactions, stoichiometry, kinetics, inlet conditions, physical prop. constants, reactor dimensions and catalyst specs) should be specified here. Constants for physical property calculations should also be updated in 'physicalData.mat' to run main_fnc. Refer to directory
/generatePhysicalData to create the data. Modify based on your system.

calculateObjective: Objective function to be minimized during the optimization is specified here. Modify based on your system.

simulations: Solves model differential equations with ode15s. Simple calculations of z,r independent variables (inlet mass flow rate, reactor dimensions) are also given here 

equations: Defines model differential equations (mass, energy and momentum balances).

calc files: Subroutines used for physical property estimations. They are called by equations.m script for ODE definitions.





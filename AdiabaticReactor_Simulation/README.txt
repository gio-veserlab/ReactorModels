This code bundle simulates a pseudohomogeneous 2D steady state tubular reactor model in adiabatic operation.Output data consists of concentration, temperature and pressure profiles across reactor, conversion and selectivity, outlet stream temperature and pressure values.

Below you can find a brief explanation of each script in the directory:

main_sim: Input data (reactions, stoichiometry, kinetics, inlet conditions, physical prop. constants, reactor dimensions and catalyst specs) should be specified here. Constants for physical property calculations should also be updated in 'physicalData.mat' to run main_fnc. Refer to directory
/generatePhysicalData to create the data. Modify based on your system.

simulations: Solves model differential equations with ode15s. Simple calculations of z,r independent variables (inlet mass flow rate, reactor dimensions) are also given here 

equations: Defines model differential equations (mass, energy and momentum balances).

calc files: Subroutines used for physical property estimations. They are called by equations.m script for ODE definitions.





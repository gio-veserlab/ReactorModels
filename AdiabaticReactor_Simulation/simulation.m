function [x, y, conversion, selectivity, delta_P, final_T, final_P, max_T, mcat, Lr, msc] ...
    = simulation(n,p,m,Rm,mw,deltaH,stoichiometry_matrix, ...
    exponent_matrix, T0, P0, tube_count, dr_h, ...
    dp, lambda_sol, n0, whsv, limiting_index, product_index, void, rho_bulk, ...
    cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, kvdip, Tbs, Vbs, mup, ...
    kr, Ea, kad, Ead)    
%% This function simulates the steady state conditions in the adiabatic reactor.

% Model parameters:

% n: Number of components in the reactor
% p: Number of reactions 
% m: Number of discrete bins for radial axis
% Rm: universal gas constant (scalar) (J/mol/K)
% mw: Molecular weights for each component (nx1) (g/mol)
% deltaH: Enthalpy of reaction (px1) (J/mol)
% stoichiometry matrix: stoichiometry of each component in each reaction (nxp) 
% exponent matrix: exponent of each component in each reaction (nxp)
% T0: inlet temperature (scalar) (C)
% P0: inlet pressure (scalar) (bar)
% tube_count: total number of tubes 
% dr_h: the diameter of the reactor (scalar) (m)
% dp: the diameter of the catalyst particles (scalar) (m)
% lambda_sol: solid bed thermal conductivity (scalar) (W m-1 K-1)
% n0: inlet mol for each component (nx1) (mol/s)
% whsv: weight hour space velocity (scalar) (mol/kg catalyst/hour)
% limiting_index: the index of the limiting reactant (1 to n)
% desiring_index: the index of the desired product (1 to n)
% void: void fraction of the catalyst (scalar) (unitless)
% rho_bulk: bulk density of catalyst (scalar) (kg/m^3)
% cpc: constants for ideal gas heat capacity for each component (nx7) (variable)
% Pcs: critical pressures for each component (nx1) (bar)
% Tcs: critical temperatures for each component (nx1) (K)
% omega: acentric factor for each component (nx1) (unitless)
% kij: binary interaction parameter for Peng-Robinson EOS (nxn) (unitless)
% Vcs: critical molar volumes (nx1) (cm3/mol)
% Zcs: critical compressibility factor (nx1) (unitless)
% muvdip: constants for ideal viscosity calculation (nx4) (variable)
% kvdip: constants for ideal thermal conductivity calculation (nx4) (variable)
% Tbs: boiling points for each component (nx1) (C)
% Vbs: molar volumes for each component at boiling points (nx1) (cm3/mol)
% mup: dipole moment for each component (nx1) (debyes)
% kr: pre-exponential factors for each reaction (1xp) (variable)
% Ea: activation energies for each reaction (1xp) (J/mol)
% kad: pre-exponential factor for adsorption kinetics (1xp) (variable)
% Ead: activation energy term for adsorption kinetics (1xp) (J/mol)

%% Calculation of the inlet mass flowrate and catalyst weight 

mIn = sum(n0 .* mw)/1000; 
% total mass flowrate, kg/s
msc = mIn / tube_count; 
% mass flowrate through single tube, kg/s
mcat = n0(limiting_index) * mw(limiting_index) / whsv * 3600 / 1000; 
% calculation of catalyst weight using whsv, kg

%% Calculation of reactor dimensions

Vcat = mcat / rho_bulk; % total catalyst volume, m3
Vr = Vcat / tube_count; % tube volume, m3
Acs = pi / 4 * dr_h^2; % tube cross sectional area, m2
Lr = Vr / Acs; % tube length, m

%%
B = 2.5*((1-void)/void)^1.11; % bed conductivity parameters for spherical
% particles, Bauer and Schlünder 

%% Input parameters to ODE solver

init = [n0/tube_count; ones(m+2,1)*T0;P0];

% n moles 11
% m+2 temperatures (m+1 grid points for m bins, + 1 average temperature for
% physical property estimation) 8
% 1 pressure 1

%% ODE integration using stiff integrator ode15s

xspan = linspace(0,Lr,51); % axial discretization for output
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[x,y] = ode15s(@(x,y) equations(x, y, m, n, p, stoichiometry_matrix, ...
    exponent_matrix, mw, Rm, deltaH, kr, Ea, kad, Ead, void, B, dr_h, ... 
    rho_bulk, msc, Acs, cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, ...
    kvdip, Tbs, Vbs, mup, dp, lambda_sol),xspan,init,opts);

%% Calculation of conversion, selectivity and outlet T-P conditions

conversion = (y(1,limiting_index) - y(end,limiting_index)) / y(1,limiting_index);
selectivity = (y(end,product_index) - y(1,product_index)) / ...
    (y(1,limiting_index) - y(end,limiting_index));

delta_P = y(1,m+n+3) - y(end,m+n+3); % pressure drop, bar
final_T = y(end,m+n+2); % outlet T, C
final_P = y(end,m+n+3); % outlet P, bar
max_T = max(y(end,m+n+2)); % maximum T, C

end

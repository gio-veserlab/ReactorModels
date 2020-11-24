function [output,  conversion, selectivity, delta_P, final_T] = ...
    calculateObjective(n,p,stoichiometry_matrix, ...
    exponent_matrix,deltaH,kr,Ea,kad,Ead,void,rho_bulk,dp,lambda_sol, ...
    n_limiting,limiting_index,ratio_index,product_index,m,Rm,mw, ...
    cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, kvdip, Tbs, Vbs, mup, ...
    limiting_ratio, T0, P0, tube_count, Lr, whsv, lambda_w, b_w, alpha_CM)

%% %% This function simulates the steady state conditions in the adiabatic 
% reactor for optimization purposes

% Model parameters:

% n: Number of components in the reactor
% p: Number of reactions
% stoichiometry matrix: stoichiometry of each component in each reaction (nxp) 
% exponent matrix: exponent of each component in each reaction (nxp)
% deltaH: Enthalpy of reaction (px1) (J/mol)
% kr: pre-exponential factors for each reaction (1xp) (variable)
% Ea: activation energies for each reaction (1xp) (J/mol)
% kad: pre-exponential factor for adsorption kinetics (1xp) (variable)
% Ead: activation energy term for adsorption kinetics (1xp) (J/mol)
% void: void fraction of the catalyst (scalar) (unitless)
% rho_bulk: bulk density of catalyst (scalar) (kg/m^3)
% dp: the diameter of the catalyst particles (scalar) (m)
% lambda_sol: solid bed thermal conductivity (scalar) (W m-1 K-1)
% n_limiting: molar flow of the limiting component (scalar) (mol/s)
% limiting_index: the index of the limiting reactant (1 to n)
% ratio_index: the index of the component given compared to limiting (1 to
% n)
% product_index: the index of the desired product (1 to n)
% m: Number of discrete bins for radial axis
% Rm: universal gas constant (scalar) (J/mol/K)
% mw: Molecular weights for each component (nx1) (g/mol)
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
% limiting ratio: the ratio of specified component to limiting component
% (scalar) (unitless)
% aHT_s: Specific heat transfer area of reactor tube (scalar)(1/m)
% lambda_w: wall thermal conductivity (scalar)(W/m/K)
% b_w: Wall thickness (scalar)(m)
% alpha_CM: Coolant overall heat transfer coefficient (scalar)(W/m2/K)
% T_cm: Coolant temperature (scalar) (C)
% T0: inlet temperature (scalar) (C)
% P0: inlet pressure (scalar) (bar)
% tube_count: total number of tubes 
% dr_h: the diameter of the reactor (scalar) (m)
% n0: inlet mol for each component (nx1) (mol/s)
% whsv: weight hour space velocity (scalar) (mol/kg catalyst/hour)
%% calculating inlet molar flows based on input parameters

% setting up inlet molar flowrates using the specified component index and
% the inlet flowrate of limiting component. 
n0 = zeros(n,1); % initialization of the molar flows 
n0(limiting_index) = n_limiting; % inlet MeOH, mol/s
n0(ratio_index) = n0(limiting_index) * limiting_ratio; % inlet H2S, mol/s

[~, ~, conversion, selectivity, delta_P, final_T, ~, ~, ~, ~, ~] ...
    = simulation_optimization(n,p,m,Rm,mw,deltaH,stoichiometry_matrix, ...
    exponent_matrix, T0, P0, tube_count, Lr, ...
    dp, lambda_sol, n0, whsv, limiting_index, product_index, void, rho_bulk, ...
    cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, kvdip, Tbs, Vbs, mup, ...
    kr, Ea, kad, Ead, lambda_w, b_w, alpha_CM);   

% the objective function is calculated based on conversion, selectivity and
% the whsv
output = -conversion*selectivity*whsv;
end
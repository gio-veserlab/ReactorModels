function [x, y, conversion, selectivity, delta_P, final_T, final_P, max_T, mcat, Lr, msc] ...
    = simulation(n,p,m,Rm,mw,deltaH,stoichiometry_matrix, ...
    exponent_matrix, T0, P0, tube_count, Lr, ...
    dp, lambda_sol, n0, whsv, limiting_index, product_index, void, rho_bulk, ...
    cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, kvdip, Tbs, Vbs, mup, ...
    kr, Ea, kad, Ead, lambda_w, b_w, alpha_CM)

%% This function simulates the steady state conditions in the polytropic 
% reactor for optimization purposes

% Model parameters:

% n: Number of components in the reactor
% p: Number of reactions 
% m: Number of discrete points in radial direction except center
% Rm: universal gas constant (scalar) (J/mol/K)
% mw: Molecular weights for each component (nx1) (g/mol)
% deltaH: Enthalpy of reaction (px1) (J/mol)
% stoichiometry matrix: stoichiometry of each component in each reaction (nxp) 
% exponent matrix: exponent of each component in each reaction (nxp)
% T0: inlet temperature (scalar) (C)
% P0: inlet pressure (scalar) (bar)
% tube_count: total number of tubes 
% Lr: length of the tube (scalar) (m)
% dp: the diameter of the catalyst particles (scalar) (m)
% lambda_sol: solid bed thermal conductivity (scalar) (W m-1 K-1)
% n0: inlet mol for each component (nx1) (mol/s)
% whsv: weight hour space velocity (scalar) (mol/kg catalyst/hour)
% limiting_index: the index of the limiting reactant (1 to n)
% desiring_index: the index of the desired product (1 to n)
% CF: thermal conductivity calculation parameter according to Bauer and 
% Schlünder (scalar) (unitless)
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
% aHT_s: Specific heat transfer area of reactor tube (scalar)(1/m)
% lambda_w: wall thermal conductivity (scalar)(W/m/K)
% b_w: Wall thickness (scalar)(m)
% alpha_CM: Coolant overall heat transfer coefficient (scalar)(W/m2/K)
% T_cm: Coolant temperature (scalar) (C)
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
Acs = Vr / Lr; % channel cross-sectional area
dr_h = (4 * Acs / pi)^0.5; % tube diameter

%% coolant and wall geometry parameters
T_cm = T0;
T_w = T_cm + 0.303; % Wall have slightly higher temperature than coolant, C
Lcr = pi * dr_h; % tube circumference, m
aHT_s = Lcr / Acs; % specific heat transfer area, m-1
Acs_w = pi * (dr_h + b_w) * b_w; % Wall cross section per channel, m2
aHT_CM = (2 * b_w + dr_h) * pi; % Specific heat transfer area of reactor tube, m-1

%%
B = 2.5*((1-void)/void)^1.11; % bed conductivity parameters for spherical
% particles, Bauer and Schlünder 

%%
% inlet parameters to the ODE solver, mole numbes multiplied by a factor
% for simplifying equations. 
init = [n0/tube_count; ones(m+2,1)*T0;P0;T_w;0]; 

% n moles 11
% m+2 temperatures (m+1 grid points for m bins, + 1 average temperature for
% physical property estimation) 8
% 1 pressure 1

xspan = linspace(0,Lr,51); 
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[x,y] = ode15s(@(x,y) equations(x, y, m, n, p, stoichiometry_matrix, ...
    exponent_matrix, mw, Rm, deltaH, kr, Ea, kad, Ead, void, B, dr_h, ... 
    rho_bulk, msc, Acs, cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, kvdip, Tbs, Vbs, mup, dp, ...
    aHT_s, lambda_w, b_w, alpha_CM, aHT_CM, T_cm, Acs_w, Lr, lambda_sol),xspan,init, opts);

%%
% The conversion and selectivity estimation based on limiting reactant and 
% product index 
conversion = (y(1,limiting_index) - y(end,limiting_index)) / y(1,limiting_index);
selectivity = (y(end,product_index) - y(1,product_index)) / ...
    (y(1,limiting_index) - y(end,limiting_index));
%%
% some final parameters from the model simulation 
delta_P = y(1,m+n+3) - y(end,m+n+3); % pressure drop, bar
final_T = y(end,n+1); % outlet temperature at the center of the tube, C
final_P = y(end,m+n+3); % outlet P, bar
max_T = max(y(end,n+1)); % maximum T, C 
end
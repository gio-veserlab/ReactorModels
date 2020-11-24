function [ ddx ] = equations(x, y, m, n, p, stoichiometry_matrix, ...
    exponent_matrix, mw, Rm, deltaH, kr, Ea, kad, Ead, void, B, dr_h, ...
    rho_bulk, msc, Acs, cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, ...
    kvdip, Tbs, Vbs, mup, dp, aHT_s, lambda_w, b_w, alpha_cm, aHT_cm, ...
    T_cm, Acs_w, Lr, lambda_sol)
%% This function calculates the derivative term for numerical integration of 
% balances for a 2D cylindrical reactor model. The model assumes the mole 
% terms does not change with radial displacement and uses method of lines to 
% define the energy balances. This function calculates n mole, m+3 energy 
% and 1 momentum balance equation. n stands for the number of components 
% in the system, m denotes the number of discrete data points in the radial
% direction to be used by method lines. For polytropic systems, there are 
% three additional derivative terms needs to be calculated. These three 
% terms are the average temperature Tavg, wall temperature Tw and the 
% derivative of the Tw which is used by the second order derivative which 
% is used to integrate wall temperature. In addition to m+n+5, there is 
% also one more equation is added for momentum balance. There are also p 
% reactions in the system.
%% PARAMETERS
% x: the location on the axial axis (m)
% y: mole, temperature and pressure values (mole: mol/s, temperature: C,
% pressure: bar)
% m: total number of discrete points on radial axis except center
% n: total number of species in the reactor
% p: total number of reactions 
% stoichiometry matrix: n by p matrix with reaction stoichiometries
% exponent matrix: n by p matrix with exponents of concentration values in
% mass action kinetics.
% mw: molecular weights of each component (nx1)(g/mol)
% Rm: universal gas constant (scalar)(J/mol/K) 
% deltaH: enthalpy of reactions (px1) (J/mol)
% kr: pre-exponential factors for each reaction (1xp) (variable)
% Ea: activation energies for each reaction (1xp) (J/mol)
% kad: pre-exponential factor for adsorption kinetics (1xp) (variable)
% Ead: activation energy term for adsorption kinetics (1xp) (J/mol)
% void: void fraction of the catalyty
% B: term used in calculation of the thermal conductivity (unitless)
% dr_h: diameter of a tube (scalar) (m)
% rho_bulk: bulk density of the catalyst (scalar) (kg/m3)
% msc: mass flow rate (scalar) (kg/s)
% Acs: cross-sectional area (scalar) (m2)
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
% dp: catalyst particle diameter (scalar) (m)
% lambda_sol: solid bed thermal conductivity (scalar) (W m-1 K-1)
% aHT_s: Specific heat transfer area of reactor tube (scalar)(1/m)
% lambda_w: wall thermal conductivity (scalar)(W/m/K)
% b_w: Wall thickness (scalar)(m)
% alpha_CM: Coolant overall heat transfer coefficient (scalar)(W/m2/K)
% aHT_CM: Specific heat transfer area of reactor (scalar)(1/m)
% T_cm: Coolant temperature (scalar) (C)
% Acs_w: Wall cross section (scalar)(m2)
% Lr: length of the reactor (scalar)(m)
%% Definitions

ddx = zeros(length(y),1); % initialization of derivatives 

n_moles = y(1:n)'; % molar flows
T_rads = y(n+1:n+m+1)'; % temperatures
Tavg = y(m+n+2); % avg temperature (for property estimation)
P = y(m+n+3); % pressure
Tw = y(m+n+4); % wall temperature
dTw = y(m+n+5); % derivative of wall temperature

%% Radial axis discretization

delta_r = dr_h / (2 * m); % finite difference spacing in r axis in m

r = linspace(0,dr_h/2,m+1); % node distance to center in m

%% Mole balance ODEs

Y = n_moles / sum(n_moles); % mole fractions at each node

[rho, ~] = calcRho(Y', P, Tavg, mw, Pcs, Tcs, omega, kij); % mixture density

v0 = msc / rho; % volumetric flow rate

C_moles = n_moles / v0; % concentration of each species

% Calculation of rate constants and adsorption constants

% the model uses Arrhenius type kinetics for reaction i:
% k_in = k_i0 * e^{-Ea/R/T} 
% for the adsorption constant the following formula is used: 
% k_ads = (1 + sigma_i=0^n exp(Ea/R/T) * C_i)^2

kin = zeros(m+1,p);
kdiv = zeros(m+1,p);

for i=1:m+1
    kin(i,:) = kr .* exp(-Ea/Rm/(T_rads(i)+273.15));
    sum_term = 0;
    for j=1:n
        sum_term = sum_term + kad(j) * exp(Ead(j)/Rm/(T_rads(i)+273.15)) * ...
            (C_moles(j) + 1e-15);
    end
    kdiv(i,:) = (ones(1,p) + ones(1,p) * sum_term).^2;
end 

% The rates of the individual reactions at different discrete points are
% calculated. 

rate = zeros(m+1,p);

for i=1:m+1 
   for j=1:p 
       rate(i,j) = kin(i,j);
       for k=1:n 
           rate(i,j) = rate(i,j) * (C_moles(k) + 1e-15)^exponent_matrix(k,j);
       end
       rate(i,j) = rate(i,j) / kdiv(i,j);
   end
end

% The average rate across r is evaluated using trapezoid method 

f_r = zeros(m+1,1); 
rateavg2 = zeros(p,1);
for i=1:p
    for j=1:m+1
        f_r(j)=4*r(j)/dr_h*rate(j,i);
    end
    rateavg2(i) = trapz(r,f_r) / max(r);
end

% Mole balance equation

for i=1:n
    ddx(i) = Acs*rho_bulk*(stoichiometry_matrix(i,:)*rateavg2 + 1e-15);
end

%% Momentum balance ODEs


u = msc / Acs / rho; % Superficial velocity per channel, m/s
% linear velocity is calculated using mass flow rate, cross-sectional area
% and density

nu_mixture = calcMixtureViscosity(Y', Tavg, muvdip, mw, P, Pcs, Tcs, ...
    omega, kij, Vcs, Zcs, mup, Vbs, Tbs); % mixture viscosity

% Ergun equation for momentum balance 

deltaP = 150 * (1-void)^2 / void^3 * nu_mixture * u / dp^2 + 1.75 * ...
    (1-void) / void^3 * rho * u^2 / dp; % Pa/m

ddx(m+n+3) = real(-deltaP / 1e5); % Pressure ODE, bar/m (1e5:unit conversion)


%% Energy balance ODEs
% Mixture thermal conductivity
lambda_mixture = calcMixtureThermalConductivity(Y', Tavg, mw, muvdip, kvdip); 

kap = lambda_sol / lambda_mixture;

% Bed thermal conductivity (Bauer and Schlünder correlation)
lambda_bed = lambda_mixture * ((1-(1-void)^0.5)+(2*(1-void)^0.5)/ ...
    (1-B/kap)*(((B*(1-1/kap))/(1-B/kap)^2)*log(kap/B)-(B-1)/(1-B/kap) - (B+1)/2));

% Mixture heat capacity 
[cp, cpres] = calcCP(Y', P, Tavg, mw, cpc, Pcs, Tcs, omega, kij);

% Mixture enthalpy
QSC = rho * cp * (Tavg + 273.15); 

% Derivative of heat capacity with respect to temperature

% Derivative of ideal gas heat capacity 
cpcs = cpc * Y';
mw_total = sum(mw .* Y');

TavgK = Tavg + 273.15;
dcpigdTavg = 2 * cpcs(2) * cpcs(3)^2 * (cpcs(3) * coth(cpcs(3) / TavgK) - TavgK) * csch(cpcs(3) / TavgK)^2 + ...
    2 * cpcs(4) * cpcs(5)^2 * (cpcs(5) * tanh(cpcs(5) / TavgK) - TavgK) * sech(cpcs(5) / TavgK)^2;
dcpigdTavg = dcpigdTavg / TavgK^4 * mw_total / 1000; % J/mol/K

% Derivative of residual heat capacity (numerical derivative is calculated
% by Taylor's approximation.)

cpress = zeros(2,1);
Ts = [Tavg Tavg+0.05]';
cpress(1) = cpres;
cpress(2) = calcCpRes(Y', P, Ts(2), Pcs, Tcs, omega, kij);
t = [(cpress(2) - cpress(1))/(Ts(2) - Ts(1))];  
dcprdTavg = t(1) * mw_total; % J/mol/K

% Total derivative of heat capacity (ideal + residual)

dcpdTavg = dcpigdTavg + dcprdTavg;

% Derivative of heat capacity with respect to pressure 

cps = zeros(2,1);
cps(1) = cp;
Ps = [P P+0.0005]';
cps(2) = calcCP(Y', Ps(2), Tavg, mw, cpc, Pcs, Tcs, omega, kij);
dcpdP = (cps(2) - cps(1))/(Ps(2) - Ps(1));

% Derivative of heat capacity with respect to molar flow of each component

dcpdn = zeros(n,1);
for i=1:n
    n_alt = n_moles';
    n_alt(i) = n_alt(i) + 1e-6;
    y_alt = n_alt / sum(n_alt);
    [cp_alt,~] = calcCP(y_alt, P, Tavg, mw, cpc, Pcs, Tcs, omega, kij);
    dcpdn(i) = (cp_alt - cp) / (n_alt(i) - n_moles(i));
end

% Derivative of density with respect to temperature 

rhos = zeros(2,1);
Ts = [Tavg Tavg+0.05]';
rhos(1) = real(rho);
rhos(2) = real(calcRho(Y',P,Ts(2),mw,Pcs,Tcs,omega,kij));
drhodTavg = (rhos(2) - rhos(1))/(Ts(2) - Ts(1));

% Derivative of density with respect to pressure 

Ps = [P P+0.0005]';
rhos(2) = real(calcRho(Y',Ps(2),Tavg,mw,Pcs,Tcs,omega,kij));
drhodP = (rhos(2) - rhos(1))/(Ps(2) - Ps(1));

% Derivative of density with respect to molar flow of each component

drhodn = zeros(n,1);
for i=1:n
    n_alt = n_moles';
    n_alt(i) = n_alt(i) + 1e-6;
    y_alt = n_alt / sum(n_alt);
    [rho_alt,~] = calcRho(y_alt, P, Tavg, mw, Pcs, Tcs, omega, kij);
    drhodn(i) = real((rho_alt - rho) / (n_alt(i) - n_moles(i)));    
end

%% Heat transfer coefficient calculation
% Reynolds and Prandtl number calculation
Re = u*dr_h*rho/nu_mixture;
Pr = nu_mixture*cp/lambda_mixture;

% Nusselt number estimation
NuxL1 = 3.66; % Nusselt number at laminar flow conditions by Lok. 
ampHT = 1.5; % Heat transfer amplification factor
if Re < 2300
    Nu = NuxL1;
elseif Re > 10000
    Nu = (((((1.8 * log10(Re) - 1.5)^(-2))/8)*Re*Pr)/(1+12.7*(((1.8 * ...
        log10(Re) - 1.5)^(-2))/8)^(0.5)*(Pr^(2./3) - 1)));
else
    Nu = (1 - (Re - 2300)/(1e4 - 2300)) * NuxL1 + ((Re - 2300)/(1e4 - ...
        2300)) * (((((5.7)^(-2))/8)*1e4*Pr)/(1+12.7*(((5.7)^(-2))/8)^(0.5)*(Pr^(2./3) - 1)));
end
alpha = ampHT * Nu * lambda_mixture / dr_h; % heat transfer coefficient

% the heat flux through wall 
q_w = alpha * aHT_s * (Tavg - Tw);

% Enthalpy balance for Tavg calculation
dummy = real((sum((rateavg2 + 1e-15) .* deltaH) * rho_bulk + q_w - ...
    ddx(m+n+3) * (-u * cp * TavgK * drhodP - u * rho * TavgK * ...
    dcpdP + QSC * u / rho * drhodP)));
for i=1:n
    dummy = dummy + ddx(i) * (-u * cp * TavgK * drhodn(i) - u * rho * ...
        TavgK * dcpdn(i) + QSC * u / rho * drhodn(i));
end
ddx(m+n+2) = dummy / (-u * cp * TavgK * drhodTavg - u * rho * cp - u * ...
    rho * TavgK * dcpdTavg + QSC * u / rho * drhodTavg);

% Wall energy balance 
% The second derivative of the wall temperature 
if x == Lr | x == 0
    ddx(m+n+4) = 0; % zero at the entrance and exit
else
    ddx(m+n+4) = y(m+n+5); % defined as another dependent variable 
end

% Energy balance ODE definition around the wall
Twa = (lambda_w * 2 * pi / log((dr_h + 2 * b_w) / dr_h) * Tw + alpha_cm * aHT_cm * T_cm) / ...
    (alpha_cm * aHT_cm + lambda_w * 2 * pi / log((dr_h + 2 * b_w) / dr_h));
ddx(m+n+5) = real((q_w * Acs / Acs_w - lambda_w * 2 * pi / log((dr_h + 2 * b_w) / dr_h) / Acs_w * (Tw - Twa)) / lambda_w);

% Energy balance equations (method of lines)
T_1 = T_rads(2); % T-1=T2 symmetry condition for center

ddx(n+1) = ((lambda_bed * ((T_1 - 2 * T_rads(1) + T_rads(2))/ delta_r^2) ...
    - sum((rateavg2 + 1e-15) .* deltaH) * rho_bulk) / (msc / Acs) - ...
    (T_rads(1) + 273.15) * (dcpdTavg *  ddx(n+m+2)+ dcpdP * ddx(n+m+3) ))/cp;
for i=2:m
   ddx(n+i) = ((lambda_bed * ((T_rads(i-1) - 2 * T_rads(i) + T_rads(i+1))/ delta_r^2 ...
       + 2 / dr_h * (T_rads(i+1) - T_rads(i-1))/delta_r) - sum((rateavg2 + 1e-15) ...
       .* deltaH) * rho_bulk) / (msc / Acs) - (T_rads(i) + 273.15) * (dcpdTavg * ...
       ddx(n+m+2)+ dcpdP * ddx(n+m+3) ))/cp;
end
ddx(n+m+1) = 0; % The temperature is constant at r=R and derivative is zero.

end


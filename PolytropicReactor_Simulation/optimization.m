% This script defines reaction system, reaction parameters, input variables,
% reactor & catalyst specifications and depending on the parameters 
% optimizes parameters.  

warning('off') % optimization can create some warning messages from ode15s
% function

% initial condition specified to begin the optimization process
ratio = 10; % initial h2s/meoh ratio
T0 = 220; % initial temp, celcius
P0 = 25; % initial pressure bar

tube_count = 1500; % initial tube count
Lr = 4.557294028; % initial tube length, m
whsv = 1; % initial space velocity

%% Definition of reactions, stoichiometry and kinetic parameters

n = 11; % Number of species
p = 11; % Number of reactions 

% 1 : ch4 
% 2 : co2 
% 3 : cs2 
% 4 : dmds
% 5 : dme 
% 6 : dms 
% 7 : h2 
% 8 : h2o 
% 9 : h2s 
% 10 : meoh
% 11 : mesh 

% Each column is a reaction and each row is a component
stoichiometry_matrix = ... 
   [0 0 0 0 0 0 0 1 1 0 0;
    0 0 0 0 0 0 0 0 1 0 0;
    0 0 0 0 0 0 0 1 0 0 0;
    0 0 0 0 0 0 0 0 0 1 -1;
    0 1 -1 -1 0 0 0 0 0 0 0;
    0 0 0 0 1 -1 1 0 0 0 0;
    0 0 0 0 0 0 0 2 2 1 -1;
    1 1 -1 0 1 0 0 0 0 0 0;
    -1 0 0 -1 0 -1 1 0 0 0 0;
    -1 -2 2 1 -1 0 0 0 -2 0 0;
    1 0 0 1 -1 2 -2 -2 0 -2 2];

% each column is a reaction and each row is a component
exponent_matrix = ...
   [0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1;
    0 0 1 0 0 0 0 0 0 0 0;
    0 0 0 0 0 1 0 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0 1;
    0 0 1 0 0 0 0 0 0 0 0;
    1 0 0 0 0 1 0 0 0 0 0;
    1 1 0 0 1 0 0 0 2 0 0;
    0 0 0 0 1 0 2 1 0 2 0];


deltaH = [-41950; -19670; 19670;  -22280; -52820; 10870; -10870; 96030; ...
    -53370; 21970; -21970]; % Enthalpy of reaction, J/mol


% Rate parameters

% pre-exponential factors
kr = [1.69E+02 9.46E-01 5.00E+00 0.00E+00 2.19E-06 2.26E+01 4.57E+00 1.90E-05 8.11E+10 8.50E-03 3.02E+00];
% activation energies
Ea = [7.70E+04 4.70E+04 6.67E+04 0.00E+00 8.00E+03 8.12E+04 7.03E+04 2.00E+04 1.68E+05 3.01E+04 8.17E+03];
% adsorption pre-exponential factors
kad = [0 0 0 0 0 0 0 5.95E-06 0 2.59E-02 4.94E-05];
% adsorption activation energies
Ead = [0 0 0 0 0 0 0 2.22E-14 0 3.11E-06 2.60E+04];
%% Input variables, reactor and catalyst specifications

dr_h = 0.0381; % tube diameter, m
rho_bulk = 801.94; % bulk catalyst density, kg/m3
dp = 0.00125; % pellet size, m
void = 0.39; % void fraction
lambda_sol = 0.2; % solid bed thermal conductivity (catalyst specific)

%% Coolant and wall material properties

lambda_w = 8; % wall thermal conductivity (W m-1 K-1)
alpha_CM = 2000; % Coolant heat transfer coefficient, W m-1 K-1
b_w = 0.013; % Wall thickness

%% Physical properties and constants

mw = [16.0428;44.0098;76.143;94.2016;46.069;62.1356;2.01588;18.0153;34.0819;32.0422;48.1088]; % molecular weight, g/mol
Rm = 8.3145; % universal gas constant

load('physicalData.mat'); 

%% Input calculations

limiting_index = 10; % limiting reactant, MeOH
product_index = 11; % main product, MeSH
ratio_index = 9; % index of reactant used to calculate from ratio (H2S)
n_limiting = 54.1857; % inlet MeOH, mol/s

%% integration variable

m = 6; % Number of discrete bins for radial axis

%% This section sets up the optimization of the defined functions with 
% given parameters and the boundaries of the parameters need to be
% optimized

% The definition of the boundaries
min_x = [9 160 15 1500 1 0.5]; % ratio temp pressure tubecount L whsv
max_x = [15.01 240.01 35.01 8000.01 8.01 6.01]; % ratio temp pressure tubecount L whsv
% Initiation of the GlobalSearch optimizer
gs = GlobalSearch('Display','iter');
% Specification of the parameters for the optimization
opts = optimoptions('fmincon','Display','iter');

% init = [12 220 25 2000 5 2];
% specifying initial parameters for optimization
init = [ratio T0 P0 tube_count Lr whsv];
% definition of the parameters need to be optimized and the objective
% function

minimize_func = @(x) calculateObjective(n,p,...
            stoichiometry_matrix, exponent_matrix,deltaH,kr,Ea,kad,Ead, ...
            void,rho_bulk,dp,lambda_sol,n_limiting,limiting_index, ...
            ratio_index,product_index,m,Rm,mw,cpc,Pcs,Tcs,omega,kij,Vcs, ...
            Zcs, muvdip, kvdip, Tbs, Vbs, mup, x(1), x(2), x(3), x(4), ...
            x(5), x(6), lambda_w, b_w, alpha_CM);
% setting the problem using the createOptimProblem function from
% optimization toolbox. 
problem = createOptimProblem('fmincon','x0',init,'objective',minimize_func,'lb',min_x,'ub',max_x,'options',opts);
% solving defined optimization problem using the global search method 
x = run(gs, problem);

% final output values are calculated using the final parameters estimated
% using global search method.
[final_values(1,1),final_values(1,2),final_values(1,3), ...
    final_values(1,4),final_values(1,5) ] = calculateObjective(n,p,...
    stoichiometry_matrix, exponent_matrix,deltaH,kr,Ea,kad,Ead, ...
    void,rho_bulk,dp,lambda_sol,n_limiting,limiting_index, ...
    ratio_index,product_index,m,Rm,mw,cpc,Pcs,Tcs,omega,kij,Vcs, ...
    Zcs, muvdip, kvdip, Tbs, Vbs, mup, ...
    x(1,1),x(1,2),x(1,3),x(1,4),x(1,5),x(1,6), ...
    lambda_w, b_w, alpha_CM); 

% printing the initial parameters and the output values
fprintf('Initial parameters\n')
fprintf('H2S/MeOH Ratio: %.2f\n',init(1));
fprintf('Inlet Temperature(C): %.2f\n',init(2));
fprintf('Inlet Pressure(bar): %.2f\n',init(3));
fprintf('Tube count: %.0f\n',init(4));
fprintf('Length of reactor(m): %.2f\n',init(5));
fprintf('WHSV: %.2f\n',init(6));
fprintf('Output parameters\n');
fprintf('H2S/MeOH Ratio: %.2f\n',x(1));
fprintf('Inlet Temperature(C): %.2f\n',x(2));
fprintf('Inlet Pressure(bar): %.2f\n',x(3));
fprintf('Tube count: %.0f\n',x(4));
fprintf('Length of reactor(m): %.2f\n',x(5));
fprintf('WHSV: %.2f\n',x(6));
fprintf('Objective Function: %.2f\n', final_values(1));
fprintf('Conversion: %.3f\n', final_values(2));
fprintf('Selectivity: %.3f\n', final_values(3));
fprintf('Pressure Drop(bar): %.2f\n', final_values(4));
fprintf('Outlet Temperature(C): %.2f\n', final_values(5));
warning('on')

% saving the output into a file
save('output.mat')
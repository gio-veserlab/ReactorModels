% This script defines reaction system, reaction parameters, input variables,
% reactor & catalyst specifications, and calls the simulation function. 

clc
%clear all
close all

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


% Rate parameters !!THE VALUES ARE CHANGED FROM DEMO TO REDACT PROPRIETARY
% INFORMATION!! 

% pre-exponential factors
kr = [1 1 1 1 1 1 1 1 1 1 1];
Ea = [1E4 1E4 1E4 1E4 1E4 1E4 1E4 1E4 1E4 1E4 1E4];
% adsorption pre-exponential factors, zero in case of no adsorption
kad = [0 0 0 0 0 0 0 0 0 0 0];
% adsorption activation energies, zero in case of no adsorption
Ead = [0 0 0 0 0 0 0 0 0 0 0];

%% Input variables, reactor and catalyst specifications

H2SMeOHratio = 10; % inlet H2S to MeOH ratio
T0 = 250; % inlet T, C
T_cm = T0; % coolant T, C
P0 = 25; % inlet P, bar
tube_count = 1500;  
dr_h = 0.0381; % tube diameter, m
whsv = 1; % residence time, h^-1
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

n0 = zeros(n,1);
n0(10) = 54.1857; % inlet MeOH, mol/s
n0(9) = n0(10) * H2SMeOHratio; % inlet H2S, mol/s

%% Model integration

m = 6; % Number of discrete bins for radial axis

[x, y, conversion, selectivity, delta_P, final_T, final_P, max_T, mcat, Lr, msc] ...
    = simulation(n,p,m,Rm,mw,deltaH,stoichiometry_matrix, ...
    exponent_matrix, T0, P0, tube_count, dr_h, ...
    dp, lambda_sol, n0, whsv, limiting_index, product_index, void, rho_bulk, ...
    cpc, Pcs, Tcs, omega, kij, Vcs, Zcs, muvdip, kvdip, Tbs, Vbs, mup, ...
    kr, Ea, kad, Ead, lambda_w, b_w, alpha_CM, T_cm);    
%% Plotting and printing, customize as needed

plot_true = 1;

%fprintf prints input & output parameters to the command window

fprintf('Input parameters\n');
fprintf('H2S/MeOH Ratio: %.2f\n',H2SMeOHratio);
fprintf('Inlet Temperature(C): %.2f\n',T0);
fprintf('Inlet Pressure(bar): %.2f\n',P0);
fprintf('Tube count: %.0f\n',tube_count);
fprintf('Length of reactor(m): %.2f\n',Lr);
fprintf('WHSV: %.2f\n',whsv);

fprintf('Output parameters\n');
fprintf('Conversion: %.3f\n', conversion);
fprintf('Selectivity: %.3f\n', selectivity);
fprintf('Productivity: %.3f\n', conversion*selectivity*whsv);


fprintf('Outlet Temperature(C): %.2f\n', final_T);
fprintf('Outlet Pressure(bar): %.2f\n', final_P);
fprintf('Pressure Drop(bar): %.2f\n', delta_P);
fprintf('Maximum Temperature(C): %.2f\n', max_T);

%figure plots for concentration, temperature and pressure profiles across z

if plot_true == 1

    close all
    set(0,'defaultAxesFontSize',18)    
    % Mole to concentration conversion
    
    sums = sum(y(:,1:11),2);
    fracs = y(:,1:11) ./ sums;
    Tavgs = mean(y(:,12:18),2);
    Tavgs = y(:,19);
    Ps = y(:,20);

    rhos = zeros(51,1);
    for i=1:51
        [rhos(i), ~] = calcRho(fracs(i,:)', Ps(i), Tavgs(i), mw, Pcs, Tcs, omega, kij);
    end

    V0s = msc ./ rhos;

    Cs = y(:,1:11) ./ V0s;

    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    % H2S
    subplot(2,3,1)
    plot(x, Cs(:,9))
    xlabel('z (m)')
    ylabel('C_{H2S} (mol/m^{3})')
    title('H2S')
    xlim([0 max(x)])
    % MESH
    subplot(2,3,2)
    plot(x, Cs(:,11))
    xlabel('z (m)')
    ylabel('C_{MeSH} (mol/m^{3})')
    title('MeSH')
    xlim([0 max(x)])
    % MEOH
    subplot(2,3,3)
    plot(x, Cs(:,10))
    xlabel('z (m)')
    ylabel('C_{MeOH} (mol/m^{3})')
    title('MeOH')
    xlim([0 max(x)])
    % DMS
    subplot(2,3,4)
    plot(x, Cs(:,6))
    xlabel('z (m)')
    ylabel('C_{DMS} (mol/m^{3})')
    title('DMS')
    xlim([0 max(x)])
    % DME
    subplot(2,3,5)
    plot(x, Cs(:,5))
    xlabel('z (m)')
    ylabel('C_{DME} (mol/m^{3})')
    title('DME')
    xlim([0 max(x)]) 
    % DMDS
    subplot(2,3,6)
    plot(x, Cs(:,4))
    xlabel('z (m)')
    ylabel('C_{DMDS} (mol/m^{3})')
    title('DMDS')
    xlim([0 max(x)])
    
    % Selectivity and conversion calculation
    
    n_meoh_out = y(end,10) ;
    n_meoh_in = n0(10)/tube_count;

    n_mesh_out = y(end,11) ;
    n_mesh_in = n0(11)/tube_count;

    n_dms_out = y(end,6) ;
    n_dms_in = n0(6)/tube_count;
    
    meoh_conversion = (n_meoh_in - n_meoh_out) / n_meoh_in;
    mesh_selectivity = (n_mesh_out - n_mesh_in) / (n_meoh_in - n_meoh_out);
    dms_selectivity = (n_dms_out - n_dms_in) / (n_meoh_in - n_meoh_out);
    
    n_meoh_out = y(:,10) ;
    n_meoh_in = n0(10)/tube_count;

    n_mesh_out = y(:,11) ;
    n_mesh_in = n0(11)/tube_count;
    
    n_dms_out = y(:,6) ;
    n_dms_in = n0(6)/tube_count;

    meoh_conversion = (n_meoh_in - n_meoh_out) / n_meoh_in;
    mesh_selectivity = (n_mesh_out - n_mesh_in) ./ (n_meoh_in - n_meoh_out(end));
    dms_selectivity = (n_dms_out - n_dms_in)*2 ./ (n_meoh_in - n_meoh_out(end));
    
    figure
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.50]);
    subplot(1,3,1)
    plot(x, y(:,12),'y')
    hold on
    plot(x, y(:,13),'r')
    plot(x, y(:,14),'g')
    plot(x, y(:,15),'k')
    plot(x, y(:,16),'m')
    plot(x, y(:,17),'c')
    plot(x, y(:,18))
    xlabel('z')
    ylabel('T(°C)')
    legend('T0','T1','T2','T3','T4','T5','T6')
    title('Temperature')
    xlim([0 max(x)])

    subplot(1,3,2)
    plot(x, y(:,20))
    xlabel('z(m)')
    ylabel('P')
    title('Pressure (bar)')
    xlim([0 max(x)])
    
    subplot(1,3,3)
    plot(x, meoh_conversion)
    hold on
    plot(x, mesh_selectivity,'r')
    plot(x, dms_selectivity,'k')
    xlabel('z(m)')
    ylabel('Conversion/Selectivity')
    title('Conversion/Selectivity')
    ax = legend('X MeOH','S MeSH','S DMS');
    set(ax,'FontSize',14)
    xlim([0 max(x)])   
    
end

function y = generatePhysicalData()

%% This function generates constants for the estimation of physical 
% properties (density, heat capacity, thermal conductivity and viscosity)
% for the species specified below. All constants are taken from Aspen
% Properties software. Please refer to physical property estimation scripts
% named "calc.." to see correlations used.

% Components:

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

Tcs = [-82.52; 31; 278.85; 341.85; 126.95; 229.89; -229.55; 374.2; 100.4; 239.43; 196.80];
Pcs = [46.17; 73.76; 79; 53.60; 53.70; 55.30; 20.47; 221.19; 90.07; 80.96; 72.3];
omega = [0.01; 0.231; 0.110697; 0.205916; 0.200221; 0.194256; 0; 0.348; 0.1; 0.557; 0.158174];
Zcs = [0.286; 0.274; 0.275; 0.264; 0.2744; 0.266; 0.305; 0.229; 0.284; 0.222; 0.268];
Vcs = [0.0986; 0.094; 0.16; 0.252; 0.17; 0.201; 0.064147; 0.0559472; 0.0985; 0.117; 0.145];

kij = zeros(11,11);
kij(10,8) = -0.0778;
kij(8,10) = kij(10,8); 
kij(9,8) = 0.04; 
kij(8,9) = kij(9,8);
kij(7,1) = 0.0156;
kij(1,7) = kij(7,1);
kij(10,2) = 0.023;
kij(2,10) = kij(10,2); 
kij(9,2) = 0.0974;
kij(2,9) = kij(9,2);
kij(8,2) = 0.12;
kij(2,8) = kij(8,2); 
kij(7,2) = -0.1622;
kij(2,7) = kij(7,2);
kij(1,2) = 0.0919;
kij(2,1) = kij(1,2); 

cpigs = readtable('cpigs.csv');
cpc = cpigs.Variables;
      
muvdip = [5.25E-07	0.59006	105.67	0;
    2.15E-06	0.46	290	0;
    5.82E-08	0.9262	44.581	0;
    3.23E-08	0.97742	0	0;
    2.68E-06	0.3975	534	0;
    5.29E-07	0.6112	302.85	0;
    1.80E-07	0.685	-0.59	140;
    1.71E-08	1.1146	0	0;
    3.93E-08	1.0134	0	0;
    3.07E-07	0.69655	205	0;
    1.64E-07	0.76706	107.97	0
    ];

kvdip = [8.40E-06	1.4268	-49.654	0;
    3.69	-0.3838	964	1.86E+06;
    0.0003467	0.7345	479	0;
    0.00022578	0.892	697	0;
    0.059975	0.2667	1018.6	1.10E+06;
    0.00023614	0.9204	638	0;
    0.002653	0.7452	12	0;
    6.20E-06	1.3973	0	0;
    1.38E-07	1.8379	-352.09	46041;
    5.80E-07	1.7862	0	0;
    2.65E-05	1.1631	29.996	32519
    ];

mup = [0; 0; 0; 1.98463; 1.3011; 1.49896; 0; 1.84972; 0.968331; 1.69983; 1.51995];
Tbs = [-161.49; -78.45; 46.225; 109.75; -24.84; 37.33; -252.76; 100; -60.35; 64.7; 5.956]; 
Vbs = [0.0379694; 0.0350189; 0.062295; 0.0980576; 0.0630445; 0.074986; 0.0285681; 0.0188311; 0.0358604; 0.0427452; 0.0542058];


save('physicalData.mat','Tcs','Pcs','omega','kij','cpc','muvdip','kvdip','Vcs','Zcs','mup','Tbs','Vbs');
y=0;
function y = generatePhysicalData()
kdata = readtable('lambda_vap.csv');
for i=1:11
   if i == 1
       pp = csapi(kdata.x___T,kdata.CH4);
   elseif i == 2
       pp = csapi(kdata.x___T,kdata.CO2);
   elseif i == 3
       pp = csapi(kdata.x___T,kdata.CS2);
   elseif i == 4
       pp = csapi(kdata.x___T,kdata.DMDS);
   elseif i == 5
       pp = csapi(kdata.x___T,kdata.DME);
   elseif i == 6
       pp = csapi(kdata.x___T,kdata.DMS);
   elseif i == 7
       pp = csapi(kdata.x___T,kdata.H2);
   elseif i == 8
       pp = csapi(kdata.x___T,kdata.H2O);
   elseif i == 9
       pp = csapi(kdata.x___T,kdata.H2S);
   elseif i == 10
       pp = csapi(kdata.x___T,kdata.MEOH);
   elseif i == 11
       pp = csapi(kdata.x___T,kdata.MESH);
   end
   ksplines{i} = pp;
end

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

    								
% 
% rhodata = readtable('rho_mass.csv');
% for i=1:11
%    if i == 1
%        pp = csapi(rhodata.x___T,rhodata.CH4);
%    elseif i == 2
%        pp = csapi(rhodata.x___T,rhodata.CO2);
%    elseif i == 3
%        pp = csapi(rhodata.x___T,rhodata.CS2);
%    elseif i == 4
%        pp = csapi(rhodata.x___T,rhodata.DMDS);
%    elseif i == 5
%        pp = csapi(rhodata.x___T,rhodata.DME);
%    elseif i == 6
%        pp = csapi(rhodata.x___T,rhodata.DMS);
%    elseif i == 7
%        pp = csapi(rhodata.x___T,rhodata.H2);
%    elseif i == 8
%        pp = csapi(rhodata.x___T,rhodata.H2O);
%    elseif i == 9
%        pp = csapi(rhodata.x___T,rhodata.H2S);
%    elseif i == 10
%        pp = csapi(rhodata.x___T,rhodata.MEOH);
%    elseif i == 11
%        pp = csapi(rhodata.x___T,rhodata.MESH);
%    end
%    rhosplines{i} = pp;
% end
% 
% cpdata = readtable('cpv_mass.csv');
% for i=1:11
%    if i == 1
%        pp = csapi(cpdata.x___T,cpdata.CH4); % *16.0428);
%    elseif i == 2
%        pp = csapi(cpdata.x___T,cpdata.CO2); % *44.0098);
%    elseif i == 3
%        pp = csapi(cpdata.x___T,cpdata.CS2); % *76.143);
%    elseif i == 4
%        pp = csapi(cpdata.x___T,cpdata.DMDS); % *94.2016);
%    elseif i == 5
%        pp = csapi(cpdata.x___T,cpdata.DME); % *46.089);
%    elseif i == 6
%        pp = csapi(cpdata.x___T,cpdata.DMS); % *62.1356);
%    elseif i == 7
%        pp = csapi(cpdata.x___T,cpdata.H2); % *2.015);
%    elseif i == 8
%        pp = csapi(cpdata.x___T,cpdata.H2O); % *18.015);
%    elseif i == 9
%        pp = csapi(cpdata.x___T,cpdata.H2S); % *34.08);
%    elseif i == 10
%        pp = csapi(cpdata.x___T,cpdata.MEOH); % *32.04);
%    elseif i == 11
%        pp = csapi(cpdata.x___T,cpdata.MESH); % *48.11);
%    end
%    cpsplines{i} = pp;
% end

save('physicalData.mat','ksplines','Tcs','Pcs','omega','kij','cpc','muvdip','kvdip','Vcs','Zcs','mup','Tbs','Vbs');
y=0;
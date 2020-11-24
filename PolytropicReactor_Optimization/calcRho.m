function [rho, z] = calcRho(y,P,T,mw,Pcs,Tcs,omega,kij)
%% calculates the mass density and molar volume based on temperature, 
% and pressure using standard peng-robinson equation of state
% y: mole fraction of each component (nx1) (unitless)
% P: pressure (scalar) (bar)
% T: temperature (scalar) (C)
% mw: molecular weights of each component (nx1) (g/mol)
% Pcs: critical pressures for each component (nx1) (bar)
% Tcs: critical temperatures for each component (nx1) (K)
% omega: acentric factor for each component (nx1) (unitless)
% kij: binary interaction parameter for Peng-Robinson EOS (nxn) (unitless)
%%
    ai = zeros(length(y),1);
    bi = zeros(length(y),1);
    % specific to PR equation the universal gas constant should be in the 
    % following units. 
    R = 8.314e-2; % L bar K-1 mol-1
    m = zeros(length(y),1);
    alpha = zeros(length(y),1);
    Tr = zeros(length(y),1);
    for i = 1:length(y)
        m(i) = 0.37464 + 1.54226 * omega(i) - 0.26992 * omega(i)^2;
        Tr(i) = (T+273.15) / (Tcs(i) + 273.15);
        alpha(i) = (1 + m(i) * (1-Tr(i)^(1/2)))^2;
        ai(i) = alpha(i) * 0.45724 * R^2 * (Tcs(i)+273.15)^2 / Pcs(i);
        bi(i) = 0.07780 * R * (Tcs(i) + 273.15) / Pcs(i);
    end
    
    % mixture parameters are calculated using mixing formulas. 
    b = sum(y .* bi);
    a = 0;
    for i=1:length(y)
        for j=1:length(y)
            a = a + y(i) * y(j) * (ai(i) * ai(j))^0.5 * (1 - kij(i,j));
        end
    end
    mw_mix = sum(mw .* y);
    % the positive solution of the Peng Robinson equation gives the molar 
    % volume in L/mol.
    res = roots([P P*b-R*(T+273.15) -3*P*b^2-2*R*(T+273.15)*b+a -a*b+R*(T+273.15)*b^2+P*b^3]);
    vm = res(1);
    % it is converted to density using mixture molecular weight.  
    rho = 1 / vm * mw_mix;
    z = P * vm / R / (T + 273.15);   
end
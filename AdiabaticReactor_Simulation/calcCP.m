function [cp cpres] = calcCP(y, P, T, mw, cpc, Pcs, Tcs, omega, kij)
%% This function calculates the heat capacity of the mixture. It is 
% calculated based on ideal gas and residual heat capacity. 
% y: mole fractions (nx1) (unitless)
% T: temperature (scalar) (C)
% P: pressure (scalar) (bar)
% cpc: constants for ideal gas heat capacity for each component (nx7) (variable)
% Pcs: critical pressures for each component (nx1) (bar)
% Tcs: critical temperatures for each component (nx1) (K)
% omega: acentric factor for each component (nx1) (unitless)
% kij: binary interaction parameter for Peng-Robinson EOS (nxn) (unitless)
%%
    n = length(y);
    
    cpigs = calcCPig(T, cpc, n); % J / kmol / K
    
    cpig = sum( y .* cpigs);
    
    cpres = calcCpRes(y, P, T, Pcs, Tcs, omega, kij);

    mw_total = sum(mw .* y);
    
    cp = (cpres + cpig) * 1000 / mw_total;
    
end 
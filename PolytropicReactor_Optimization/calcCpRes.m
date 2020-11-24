function cpres = calcCpRes(y, P, T, Pcs, Tcs, omega, kij)
%% This function calculates the residual heat capacity by taking numerical 
% derivative of residual enthalpy calculated from Peng-Robinson equation
% of state with respect to temperature change.
% y: mole fractions (nx1) (unitless)
% T: temperature (scalar) (C)
% P: pressure (scalar) (bar)
% cpc: constants for ideal gas heat capacity for each component (nx7) (variable)
% Pcs: critical pressures for each component (nx1) (bar)
% Tcs: critical temperatures for each component (nx1) (K)
% omega: acentric factor for each component (nx1) (unitless)
% kij: binary interaction parameter for Peng-Robinson EOS (nxn) (unitless)
%%
    Hrs = zeros(2,1);
    Ts = [T T+0.05]';
    for i=1:length(Ts)
        Hrs(i) = real(calcHr(y, P, Ts(i), Pcs, Tcs, omega, kij));
    end
    p = [(Hrs(2) - Hrs(1))/(Ts(2) - Ts(1))];
    
    cpres = p(1);
    
end
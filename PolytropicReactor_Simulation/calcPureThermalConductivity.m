function k_pure_zero = calcPureThermalConductivity(T ,kvdip)
%% calculates the pure component thermal conductivity for all the components
% in N * s / m^2
% T: temperature (scalar) (C)
% kvdip: constants for pure viscosity calculation by DIPPR methods
%%
    n = size(kvdip,1);
    k_pure_zero = zeros(n,1);
    
    for i=1:n
        k_pure_zero(i,1) = kvdip(i,1) * (T+273.15)^kvdip(i,2) / (1 + kvdip(i,3) / (T+273.15) + kvdip(i,4) / (T + 273.15)^2 );  % N * s / m^2
    end
end
function cpigs = calcCPig(T, cpc, n)
% This function calculates the ideal gas heat capacity using DIPPR formula.
% T: temperature (scalar) (C)
% cpc: constants for ideal gas heat capacity for each component (nx7) (variable)
% n: number of components
    cpigs = zeros(n,1);
    for i=1:n
        cpigs(i,1) = cpc(1,i) + cpc(2,i) * ((cpc(3,i) / (T + 273.15)) / sinh(cpc(3,i) / (T + 273.15)))^2 + ...
            + cpc(4,i) * ((cpc(5,i) / (T + 273.15)) / cosh(cpc(5,i) / (T + 273.15)))^2;
        cpigs(i,1) = cpigs(i,1) / 1e3;  % J / mol / K
    end
end
function k_mixture = calcMixtureThermalConductivity(y, T, mw, muvdip, kvdip)
%% This function calculates the thermal conductivity of mixture at given 
% temperature and composition. 
% y: mole fractions (nx1) (unitless)
% T: temperature (scalar) (C)
% mw: molecular weights (nx1) (g/mol)
% muvdip: constants for ideal viscosity calculation (nx4) (variable)
% kvdip: constants for ideal thermal conductivity calculation (nx4) (variable)

    viscosity = calcPureViscosity(T, muvdip);
    n = length(y);
    Aij = zeros(n,n);
    % mixing parameters for thermal conductivity is calculated using the
    % pure component viscosities at low temperature and pressures. 
    for i = 1:n
        for j = 1:n
            Aij(i,j) = (1 + (viscosity(i,1) / viscosity(j,1))^0.5 * ...
                (mw(j) / mw(i))^0.25)^2 / (8 * (1 + mw(i) / mw(j)))^0.5;
        end
    end
    
    k_pures = calcPureThermalConductivity(T ,kvdip);
    % mixture thermal conductivity is calculated using pure component
    % thermal conductivity and binary interaction
    % parameters evaluated by pure component viscosity
    k_mixture = 0;
    for i = 1:n
        summation = 0;
        for j = 1:n
            summation = summation + y(j) * Aij(i,j);
        end
        k_mixture = k_mixture + y(i) * k_pures(i) / summation;
    end
end

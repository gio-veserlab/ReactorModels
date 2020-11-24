function nu = calcPureViscosity(T,muvdip)
%% calculates the pure pure component viscosity for all the components
% in N * s / m^2
% T: temperature (scalar) (C)
% muvdip: constants for pure viscosity calculation by DIPPR methods
%%
    n = size(muvdip,1);
    nu = zeros(n,1);
    for i=1:n
        nu(i,1) = muvdip(i,1) * (T + 273.15)^muvdip(i,2) / (1 + muvdip(i,3) / (T + 273.15) + muvdip(i,4) / (T + 273.15)^2 );  % N * s / m^2
    end
end
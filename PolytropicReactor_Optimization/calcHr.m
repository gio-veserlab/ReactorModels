function HR = calcHr(y, P, T, Pcs, Tcs, omega, kij)
%% This function calculates the residual enthalpy of a given mixture using
% Peng-Robinson equation of state. 
% y: mole fractions (nx1) (unitless)
% T: temperature (scalar) (C)
% P: pressure (scalar) (bar)
% cpc: constants for ideal gas heat capacity for each component (nx7) (variable)
% Pcs: critical pressures for each component (nx1) (bar)
% Tcs: critical temperatures for each component (nx1) (K)
% omega: acentric factor for each component (nx1) (unitless)
% kij: binary interaction parameter for Peng-Robinson EOS (nxn) (unitless)
%%

    e  = 1 - sqrt(2);
    s  = 1 + sqrt(2);
    ai = zeros(length(y),1);
    bi = zeros(length(y),1);
    m = zeros(length(y),1);
    alpha = zeros(length(y),1);
    R = 8.314e-2; % L bar K-1 mol-1
    Tr = zeros(length(y),1);
    for i = 1:length(y)
        m(i) = 0.37464 + 1.54226 * omega(i) - 0.26992 * omega(i)^2; %% Unitless
        Tr(i) = (T+273.15) / (Tcs(i) + 273.15); % Unitless 
        alpha(i) = (1 + m(i) * (1-Tr(i)^(1/2)))^2; %% Unitless
        ai(i) = alpha(i) * 0.45724 * R^2 * (Tcs(i)+273.15)^2 / Pcs(i);  % L^2 bar / mol^2
        bi(i) = 0.07780 * R * (Tcs(i) + 273.15) / Pcs(i); % L / mol
    end
    
    b = sum(y .* bi); % L / mol
    a = 0; % L^2 bar / mol^2
    for i=1:length(y)
        for j=1:length(y)
            a = a + y(i) * y(j) * (ai(i) * ai(j))^0.5 * (1 - kij(i,j));
        end
    end
    
    Q = ((ai*ai').^0.5).*(1 - kij);
    dQdT =  0.45724*(R^2)*(kij - 1).*((Tcs+273.15)*(Tcs+273.15)')./((Pcs*Pcs').^0.5).*(1/(2*(T+273.15)^0.5)).*...
    ((m./((Tcs+273.15).^0.5))*(alpha.^0.5)' + (alpha.^0.5)*(m./((Tcs+273.15).^0.5))');

    dadT = y'*(Q - (T + 273.15)*dQdT)*y;
    res = roots([P P*b-R*(T+273.15) -3*P*b^2-2*R*(T+273.15)*b+a -a*b+R*(T+273.15)*b^2+P*b^3]);
    vm = res(1);
    z = P*vm/(R*(T+273.15));
    HR  = (R*(T+273.15)*(z - 1) - dadT/(2*(s - 1)*b)*log((vm + s*b)/(vm + e*b))) * 100; % J/mol
end
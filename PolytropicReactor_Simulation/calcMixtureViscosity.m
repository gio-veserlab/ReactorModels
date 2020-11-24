function nu_mixture = calcMixtureViscosity(y, T, muvdip, mw, P, Pcs, Tcs, omega, kij, Vcs, Zcs, mup, Vbs, Tbs)
%% This function calculates the viscosity of mixture at given 
% temperature and composition. 
% y: mole fractions (nx1) (unitless)
% T: temperature (scalar) (C)
% muvdip: constants for ideal viscosity calculation (nx4) (variable)
% mw: molecular weights (nx1) (g/mol)
% P: pressure (scalar) (bar)
% Pcs: critical pressures for each component (nx1) (bar)
% Tcs: critical temperatures for each component (nx1) (K)
% omega: acentric factor for each component (nx1) (unitless)
% kij: binary interaction parameter for Peng-Robinson EOS (nxn) (unitless)
% Vcs: critical molar volumes (nx1) (cm3/mol)
% Zcs: critical compressibility factor (nx1) (unitless)
% mup: dipole moment for each component (nx1) (debyes)
% Tbs: boiling points for each component (nx1) (C)
% Vbs: molar volumes for each component at boiling points (nx1) (cm3/mol)
%%
    nu_pure = calcPureViscosity(T, muvdip);
    n = length(y);
    % boiling values and dipole moment is used to estimate normalized
    % temperature (T_star).
    epsilon_k = zeros(n,1);
    T_star = zeros(n,1);
    del = zeros(n,1);
    for i=1:n
        del(i,1) = 2 * mup(i,1)^2 / Vbs(i,1) / (Tbs(i,1) + 273.15);
        epsilon_k(i,1) = 1.15 * (1 + 0.85*del(i,1)^2) * (Tbs(i,1) + 273.15);
        T_star(i,1) = (T + 273.15) / epsilon_k(i,1);
    end
    
    % various binary parameters are calculated. The final one is the phi
    % parameter. 
    S_ij = zeros(n,n);
    for i=1:n
        for j=1:n
            S_ij(i,j) = (1 + (T_star(i,1)*T_star(j,1))^.5 + del(i,1) * del(j,1) / 4) / ((1 + T_star(i,1) + del(i,1)^2 / 4)^.5 * (1 + T_star(j,1) + del(j,1)^2 / 4)^.5);
        end
    end
    
    M_ij = zeros(n,n);
    m_ij = zeros(n,n);
    A_ij = zeros(n,n);
    for i=1:n
        for j=1:n
            M_ij(i,j) = mw(i) / mw(j);
            m_ij(i,j) = (4 / (1+M_ij(i,j)^-1) / (1 + M_ij(i,j)))^(1./4);
            A_ij(i,j) = m_ij(i,j) * M_ij(i,j)^(-.5) * (1 + (M_ij(i,j) - M_ij(i,j)^.45) / (2 * (1 + M_ij(i,j))+ ((1 + M_ij(i,j)^.45) * m_ij(i,j))/(1 + m_ij(i,j))));
        end
    end        
    
    phi_ij=zeros(n,n);
    for i=1:n
        for j=1:n
            phi_ij(i,j) = (nu_pure(i) / nu_pure(j))^0.5 * S_ij(i,j) * A_ij(i,j);
        end
    end
        
    % After binary interaction parameters estimated, the mixture viscosity
    % at low pressure is calculated. 
    nu_lp_mixture = 0;
    for i=1:n
        summation = 0;
        for j=1:n
            if i ~= j
                summation = summation + y(j) * phi_ij(i,j);
            end
        end
        nu_lp_mixture = nu_lp_mixture + y(i) * nu_pure(i) / (y(i) + summation);
    end

    % Temperature correction requires the reduced molar volume of the
    % mixture. It is calculated using Peng-Robinson equation of state. 
    R = 8.314;
    [~, z] = calcRho(y,P,T,mw,Pcs,Tcs,omega,kij);
    v = z * R * 1e-2 * (T + 273.15) / P;
    Vc = sum(y .* Vcs);
    rho_r = Vc / v;
    Zc = sum(y .* Zcs);
    mw_mixture = sum(y .* mw);
    Tc = sum(y .* (Tcs + 273.15));
    Pc = Zc * R * 1e-2 * Tc / Vc;
    ksi = Tc^(1.0/6) / mw_mixture^(1./2) / Pc^(2./3);
    
    % Depending on the value of reduced density, the pressure correction is
    % performed on the low pressure mixture viscosity. 
    if rho_r < 0.1
        nu_mixture = nu_lp_mixture * 1e7 + 1.656 * rho_r^1.1111 / ksi;
    elseif rho_r > 0.1 && rho_r < 0.9
        nu_mixture = nu_lp_mixture * 1e7 + (0.0607 * (9.045 * rho_r + 0.63)^1.739) / ksi;
    elseif rho_r > 0.9 && rho_r < 2.6
        if rho_r < 2.2
            delta = 0;
        else
            delta = (4.75e-4)*(rho_r^3 - 10.65)^2;
        end
        nu_mixture = nu_lp_mixture * 1e7 + exp(4 - exp(0.6439 - 0.1005 * rho_r - delta)) / ksi;
    else
        nu_mixture = nu_lp_mixture * 1e7;
    end
    nu_mixture = nu_mixture * 1e-7;
end

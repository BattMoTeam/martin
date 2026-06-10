% Free energy f(c) = Nv [omega*(1-2c) + kT*ln(c/(1-c))]
% normalised for now (without the Nv prefactor) 

function f = computePhaseFieldEnergy(c, omega, kT)

    
    % c = max(c, 1e-10);
    % c = min(c, 1 - 1e-10);

    % epsi = 1e-12;
    % f = omega .* (1 - 2 .* c) + kT .* log((c + epsi) ./ (1 - (c - epsi)));

    f = omega .* (1 - 2 .* c) + kT .* log(c ./ (1 - c));
    
end

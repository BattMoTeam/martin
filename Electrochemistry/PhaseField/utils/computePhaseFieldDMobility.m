% Derivative of mobility M'(c) = 1-2c

function dmob = computePhaseFieldDMobility(c, L0)
% L0 is kinetic coefficient, see Han
    
    dmob = L0*(1 - 2 .* c);

end

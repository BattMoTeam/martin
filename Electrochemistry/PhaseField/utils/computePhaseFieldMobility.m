function mob = computePhaseFieldMobility(c, L0)
% L0 is kinetic coefficient, see Han
    mob = L0*c.*(1 - c);
    
end


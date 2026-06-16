classdef PhaseFieldInputParams < InputParams 
%
% Inputs parameters for the Cahn-Hilliard phase field model
%
    properties

        %% Discretization parameters

        % Number of discretization intervals for the spectral method [-]
        N

        %% Physical functions 

        % Mobility function M(c)
        mobility
        % Derivative of mobility function M'(c)
        dMobility
        % Free energy F(c) 
        energy

        % Kinetic coefficient in m^2/J
        kineticCoefficient  % denoted L0 in ref1

        % Particle radius in m
        particleRadius

        % Interface width parameter
        epsilon

        omega % constant in the expression of the energy
              % is related to the nearest-neighbor interaction strength between lithium ions within the host
        kT    % appears in the expression of the energy
              % value normalised for now

        saturationConcentration % the saturation concentration of the guest molecule in the host material (symbol: cmax)
        
        % Advanced parameters
        np % Number of particles (will be set when initialized from above)
        volumeFraction % 
        
    end
    
    methods
        
        function inputparams = PhaseFieldInputParams(jsonstruct)
            inputparams = inputparams@InputParams(jsonstruct);
        end
        
    end
    
    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}

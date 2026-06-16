classdef PhaseField < BaseModel
% @article{ref1,
%   title   = {Electrochemical modeling of intercalation processes with phase field models},
%   journal = {Electrochimica Acta},
%   author  = {Han, B.C. and Van der Ven, A. and Morgan, D. and Ceder, G.},
%   year    = 2004,
% }
    
    properties

        %% Input parameters

        % Standard parameters

        N         % discretization parameter (size of coefC and coefW is N + 1)
        mobility  % mobility function
        dMobility % derivative of mobility function
        energy    % Free energy function

        kineticCoefficient  % denoted L0 in ref1 (used as scaling coefficient for mobility)
        particleRadius  
        
        epsilon % interface width parameter

        omega % constant in the expression of the energy
              % is related to the nearest-neighbor interaction strength between lithium ions within the host
        kT    % appears in the expression of the energy
              % value normalised for now
        
        saturationConcentration % the saturation concentration of the guest molecule in the host material (symbol: cmax)
        
        % Advanced parameters
        np % Number of particles (will be set when initialized from above)
        volumeFraction
        
        %% Helper structures
        
        A   % Collocation mapping matrices
        dA  % Collocation mapping matrices, first derivative
        ddA % Collocation mapping matrices, second derivative
        
        mobilityFunc  % mobility function
        dMobilityFunc % derivative of mobility function
        energyFunc    % Free energy function

        wScaling % scaling value for w (also called w0 in code)
        
    end
    
    methods
        
        function model = PhaseField(inputparams)
        %
        % ``inputparams`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@BaseModel();

            fdnames = {'N'                 , ...
                       'np'                , ...
                       'volumeFraction'    , ...
                       'epsilon'           , ...
                       'mobility'          , ...
                       'dMobility'         , ...
                       'energy'            , ...
                       'kineticCoefficient', ...
                       'particleRadius'    , ...
                       'omega'             , ...
                       'kT'                , ...
                       'saturationConcentration'};
            
            model = dispatchParams(model, inputparams, fdnames);

            model.operators = model.setupOperators();
            
            % when more parameters are used 
            % 
            % func  = setupFunction(model.mobility)
            % model.mobilityFunc = @(c) func(c, model.p1, model.p2);

            func  = setupFunction(model.mobility);
            model.mobilityFunc  = @(c) func(c, model.kineticCoefficient);
            
            func  = setupFunction(model.dMobility);
            model.dMobilityFunc = @(c) func(c, model.kineticCoefficient);
            
            func                = setupFunction(model.energy);
            model.energyFunc    = @(c, T) func(c, T, model.omega);
            
            model = model.setupSpectralModel();

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};
            % Spectral decomposition coefficents for c
            varnames{end + 1} = 'coefC';
            % Spectral decomposition coefficents for w
            varnames{end + 1} = 'coefW';
            % collocation values for c
            varnames{end + 1} = 'c';
            % collocation values for w
            varnames{end + 1} = 'w';
            % collocation values for dC
            varnames{end + 1} = 'dC';
            % collocation values for dW
            varnames{end + 1} = 'dW';
            % collocation values for ddC
            varnames{end + 1} = 'ddC';
            % collocation values for ddW
            varnames{end + 1} = 'ddW';
            % equation for c
            varnames{end + 1} = 'eqC';
            % equation for w
            varnames{end + 1} = 'eqW';
            % accumulation term for c
            varnames{end + 1} = 'massAccumC';
            % boundary flux at the right boundary (outward is positive)
            varnames{end + 1} = 'bdOutFlux';
            % temperature
            varnames{end + 1} = 'T';
            % average concentration mol/m^3
            varnames{end + 1} = 'cAverage';
            % surface concentration mol/m^3
            varnames{end + 1} = 'cSurface';
            
            model = model.registerVarNames(varnames);

            if model.isRootSimulationModel

                model = model.registerVarName('time');
                model = model.setAsStaticVarName('time');
                model = model.setAsExtraVarName('cSurface');
                
            else
                
                varnames = {};
                varnames{end + 1} = 'Rvol';
                model = model.registerVarNames(varnames);

            end
            
            model = model.setAsExtraVarName('cAverage');

            if model.isRootSimulationModel
                
                inputnames = {'time'};
                fn = @ProtonicMembrane.updateBdFluxFromTime;
                fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
                model = model.registerPropFunction({'bdOutFlux', fn, inputnames});
                model = model.registerPropFunction({'T', fn, inputnames});
                
            else
                fn = @PhaseField.updateBdFlux;
                model = model.registerPropFunction({'bdOutFlux', fn, {'Rvol'}});
            end

            fn = @PhaseField.updateAverageConcentration;
            model = model.registerPropFunction({'cAverage', fn, {'c'}});
            
            fn = @PhaseField.updateCsurface;
            model = model.registerPropFunction({'cSurface', fn, {'c'}});
            
            fn = @PhaseField.updateC;
            model = model.registerPropFunction({'c', fn, {'coefC'}});
            
            fn = @PhaseField.updateW;
            model = model.registerPropFunction({'w', fn, {'coefW'}});

            fn = @PhaseField.updateDC;
            model = model.registerPropFunction({'dC', fn, {'coefC'}});
            
            fn = @PhaseField.updateDW;
            model = model.registerPropFunction({'dW', fn, {'coefW'}});

            fn = @PhaseField.updateDDC;
            model = model.registerPropFunction({'ddC', fn, {'coefC'}});
            
            fn = @PhaseField.updateDDW;
            model = model.registerPropFunction({'ddW', fn, {'coefW'}});
            
            fn = @PhaseField.updateEqC;
            model = model.registerPropFunction({'eqC', fn, {'massAccumC', 'c', 'dC', 'dW', 'ddW'}});

            fn = @PhaseField.updateEqW;
            model = model.registerPropFunction({'eqW', fn, {'w', 'ddC', 'c', 'T', 'bdOutFlux'}});

            fn = @PhaseField.updateMassAccumC;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccumC', fn, {'c'}});
            
        end
        
        function operators = setupOperators(model)

            np  = model.np;
            nPt = model.N + 1;

            indInnerBc = (1 : np)';
            indInnerBc = 1 + nPt*(indInnerBc - 1);
            
            indBc = nPt*(1 : np)';

            operators = struct('indBc', indBc, ...
                               'indInnerBc', indInnerBc);

        end

        function model = setupForSimulation(model)
            
            model = model.equipModelForComputation();

            % model = model.setupScalings([]);
            
        end

        function forces = getValidDrivingForces(model)
            % needed by MRST
            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];

        end

        function newstate = addVariablesAfterConvergence(model, newstate, state)

            newstate.c = state.c;
            
        end

        function [c1, c2] = getEquilibriumValues(model, T)

            func = @(c) model.energyFunc(c, T);

            % compute concentrations at the 'dips' of the double well curve
            c0 = 0.2;
            c1 = fzero(func, c0);
            
            c0 = 0.8;
            c2 = fzero(func, c0);
            
        end

        function x = chebyshevNodes(model, N)
            % returns the N+1 Chebyshev-Lobatto nodes on [-1, 1] : x_j = cos(pi*j/N), j=0..N
            % sorted in ascending order from -1 to 1
            
            jIdx = (0 : N)';
            
            x = flipud(cos(pi * jIdx / N));

        end

        function model = setupSpectralModel(model)
            % assign collocation mappings matrices (A, dA, ddA)
            %
            % spectral decomposition : c(x,t) = sum_{i=0}^{N} coefC_i(t) * T_i(x)
            % with T_n the nth Chebyshev polynomial of the first kind
            %
            % A(j,i)   = T_{i-1}(x_{j-1})    ->  c    = A   * coefC 
            % dA(j,i)  = T'_{i-1}(x_{j-1})   ->  dC   = dA  * coefC
            % ddA(j,i) = T''_{i-1}(x_{j-1})  ->  ddC  = ddA * coefC

            x = model.chebyshevNodes(model.N);      % column vector (N+1)x1

            nPt = model.N + 1;

            chebA   = zeros(nPt, nPt);
            chebDA  = zeros(nPt, nPt);
            chebDDA = zeros(nPt, nPt);

            % init: degrees 0 and 1 (columns 1 and 2)
            chebA(:, 1)   = 1;              % T_0(x)   = 1
            chebA(:, 2)   = x;              % T_1(x)   = x

            chebDA(:, 1)  = 0;              % T'_0(x)  = 0
            chebDA(:, 2)  = 1;              % T'_1(x)  = 1

            chebDDA(:, 1) = 0;              % T''_0(x) = 0
            chebDDA(:, 2) = 0;              % T''_1(x) = 0

            % recurrence for degrees 2..N (columns 3..N+1)
            for iIdx = 3:nPt
                chebA(:, iIdx)   = 2 .* x .* chebA(:, iIdx-1)   - chebA(:, iIdx-2);
                chebDA(:, iIdx)  = 2 .* chebA(:, iIdx-1)  + 2 .* x .* chebDA(:, iIdx-1)  - chebDA(:, iIdx-2);
                chebDDA(:, iIdx) = 4 .* chebDA(:, iIdx-1) + 2 .* x .* chebDDA(:, iIdx-1) - chebDDA(:, iIdx-2);
            end

            np = model.np;

            chebA   = repmat({chebA}, 1, np);
            chebDA  = repmat({chebDA}, 1, np);
            chebDDA = repmat({chebDDA}, 1, np);
            
            % assign to model properties
            model.A   = blkdiag(chebA{:});
            model.dA  = blkdiag(chebDA{:});
            model.ddA = blkdiag(chebDDA{:});

        end

        function initstate = setupInitialState(model)

            %
            Tinit = model.kT/PhysicalConstants.kb;
            
            % initialize value of concentration
            nPt = model.N + 1;

            np = model.np;
            % 1. ------------------------------------------------
            % with a random perturbation aroud the mean value 0.5
            for ip = 1 : np
                rng(ip); % to keep the same random perturbation
                mean = 0.50;
                cInit{ip} = mean + 0.02 * randn(nPt, 1);
            end
            cInit = vertcat(cInit{:});

            % compute initial coefficients of c by inverting matrix A
            coefCInit = model.A \ cInit;

            % % 2. ---------------------------------------------------
            % % smooth perturbation : we define the coefficients first
            % rng(10); % to keep the same random perturbation
            % mean = 0.3;
            % coefCInit    = zeros(nPt, 1);
            % coefCInit(1) = mean;                   % T_0 : mean concentration
            % nModes       = 10;                     % only first modes
            % coefCInit(2 : nModes + 1) = 0.01 * randn(nModes, 1);
            % % recompute cInit from the spectral coefficients
            % cInit = model.A * coefCInit;


            % compute ddC to compute w0
            initstate.coefC = coefCInit;
            initstate = model.updateDDC(initstate);
            ddCInit = initstate.ddC;
            
            wInit = model.energyFunc(cInit, Tinit) - model.epsilon^2 .* ddCInit;

            % compute initial coefficients of w 
            coefWInit = model.A \ wInit;
            % coefWInit = zeros(nPt, 1);
            
            % initialize primary variables
            initstate.coefC = coefCInit;
            initstate.coefW = coefWInit;

            % initialize other variables
            initstate.c   = cInit;
            initstate.w   = wInit;

        end
        
        function state = updateBdFluxFromTime(model, state, drivingForces)
            % update the flux boundary condition at x=1 

            time = state.time;
            state.bdOutFlux = drivingForces.src(time);
            state.T = model.kT/PhysicalConstants.kb;
            
        end

        function state = updateCsurface(model, state)
        % update cSurface
        % NOTE : cSurface is mol/m^3 while c is without unit

            op = model.operators;
            
            state.cSurface = cmax*state.c(op.indBc);
            
        end

        
        function state = updateBdFlux(model, state)

            vf   = model.volumeFraction;
            rp   = model.particleRadius;
            cmax = model.saturationConcentration;
            
            Rvol = state.Rvol;
            
            state.bdOutFlux = (rp/(vf*cmax))*Rvol;
            
        end
        
        function state = updateC(model, state)

            state.c = model.A * state.coefC;
            
        end

        function state = updateW(model, state)

            state.w = model.A * state.coefW;

        end
        
        function state = updateDC(model, state)

            rp = model.particleRadius;
            
            state.dC = (2/rp) * model.dA * state.coefC;

        end
        
        function state = updateDW(model, state)

            rp = model.particleRadius;

            state.dW = (2/rp) * model.dA * state.coefW;

        end
        
        function state = updateDDC(model, state)

            rp = model.particleRadius;
            
            state.ddC = (2/rp)^2 * model.ddA * state.coefC;

        end
        
        function state = updateDDW(model, state)

            rp = model.particleRadius;
            
            state.ddW = (2/rp)^2 * model.ddA * state.coefW;
            
        end


        function state = updateMassAccumC(model, state, state0, dt)

            state.massAccumC = (1/dt) .* (state.c - state0.c);

        end
        
        function state = updateEqC(model, state)
            % Residual form of :
            % dc/dt - dMobility(c) * dc/dx * dw/dx - Mobility(c) * d2w/dx2 = 0

            
            massAccumC = state.massAccumC;
            c          = state.c;
            dC         = state.dC;
            dW         = state.dW;
            ddW        = state.ddW;

            mobility  = model.mobilityFunc(c);
            dMobility = model.dMobilityFunc(c);
            w0        = model.wScaling;
            rp        = model.particleRadius;
            m0        = model.kineticCoefficient; % scaling for the mobility
            op        = model.operators;

            % scaling for this equation is (w0*m0)/(rp^2)
            eqC = massAccumC - dMobility .* dC .* dW - mobility .* ddW;

            % Neumman boundary condition : dc/dn = 0
            % -> flat concentration profile at the boundary
            % 1D case : dc/dx = 0
            scaling = ((w0*m0)/(rp^2))/(1/rp);
            eqC(op.indInnerBc) = scaling*dC(op.indInnerBc);
            eqC(op.indBc)      = scaling*dC(op.indBc);
            
            state.eqC = eqC;

        end
        
        function state = updateEqW(model, state)
            % Residual form of : w + d2c/dx2 - F(c) = 0

            w0       = model.wScaling;
            m0       = model.kineticCoefficient; % scaling for mobility
            rp       = model.particleRadius;
            op       = model.operators;
            epsi     = model.epsilon;

            c         = state.c;
            T         = state.T;
            ddC       = state.ddC;
            w         = state.w;
            dW        = state.dW;
            bdOutFlux = state.bdOutFlux;
            
            F        = model.energyFunc(c, T);
            mobility = model.mobilityFunc(c);
            
            % scaling of equation is w0
            eqW = w + epsi^2 .* ddC - F;

            % Neumann boundary condition 
            % flux is given by J = -M(c) * grad(w) = -M(c) * dw/dx in 1D
            % we want the residual to be M(c) * dw/dx - J = 0
            scaling = w0/(m0*(w0/rp));
            eqW(op.indInnerBc) = scaling*(mobility(op.indInnerBc) .* dW(op.indInnerBc));
            eqW(op.indBc)      = scaling*(mobility(op.indBc) .* dW(op.indBc) + bdOutFlux);

            state.eqW = eqW;
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            % cmin = model.cmin;
            % state.(elyte).c = max(cmin, state.(elyte).c);

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

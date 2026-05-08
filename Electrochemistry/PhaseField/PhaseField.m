classdef PhaseField < BaseModel
    
    properties

        %% Input parameters

        % Standard parameters

        N         % discretization parameter
        mobility  % mobility function
        dMobility % derivative of mobility function
        energy    % Free energy function

        epsilon % interface width parameter
        
        % boundaryConditionType % 'neumann' or 'dirichlet'
        % boundaryValue         % boundary values [left, right]
        
        
        %% Helper structures
        
        A   % Collocation mapping matrices
        dA  % Collocation mapping matrices, first derivative
        ddA % Collocation mapping matrices, second derivative
        
        mobilityFunc  % mobility function
        dMobilityFunc % derivative of mobility function
        energyFunc    % Free energy function

    end
    
    methods
        
        function model = PhaseField(inputparams)
        %
        % ``inputparams`` is instance of :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
        %
            model = model@BaseModel();

            fdnames = {'N'      , ...
                       'epsilon', ...
                       'mobility'};
            
            % fdnames = {'N'                    , ...
            %            'epsilon'              , ...
            %            'mobilityFunc'         , ...
            %            'dMobilityFunc'        , ...
            %            'energyFunc'           , ...
            %            'boundaryConditionType', ...
            %            'boundaryValue'};

            model = dispatchParams(model, inputparams, fdnames);


            % func  = setupFunction(model.mobility)
            % model.mobilityFunc = @(c) func(c, model.p1, model.p2);

            model.mobilityFunc = setupFunction(model.mobility);
            
            % model = model.setupSpectralModel();

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
            

            model = model.registerVarNames(varnames);

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
            model = model.registerPropFunction({'eqW', fn, {'w', 'ddC', 'c'}});

            fn = @PhaseField.updateMassAccumC;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccumC', fn, {'c'}});
            
        end

        function model = setupForSimulation(model)
            
            model = model.equipModelForComputation();
            % model = model.setupScalings([]);
            % 
        end


        function x = chebyshevNodes(model, N)
            % returns the N+1 Chebyshev-Lobatto nodes on [-1, 1] : x_j = cos(pi*j/N), j=0..N
            % sorted in ascending order from -1 to 1
            
            jIdx = (0:N)';
            
            % size(jIdx)
            % size(N)

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

            numPoints = model.N + 1;

            chebA   = zeros(numPoints, numPoints);
            chebDA  = zeros(numPoints, numPoints);
            chebDDA = zeros(numPoints, numPoints);

            % init: degrees 0 and 1 (columns 1 and 2)
            chebA(:, 1)   = 1;              % T_0(x)   = 1
            chebA(:, 2)   = x;              % T_1(x)   = x

            chebDA(:, 1)  = 0;              % T'_0(x)  = 0
            chebDA(:, 2)  = 1;              % T'_1(x)  = 1

            chebDDA(:, 1) = 0;              % T''_0(x) = 0
            chebDDA(:, 2) = 0;              % T''_1(x) = 0

            % recurrence for degrees 2..N (columns 3..N+1)
            for iIdx = 3:numPoints
                chebA(:, iIdx)   = 2 .* x .* chebA(:, iIdx-1)   - chebA(:, iIdx-2);
                chebDA(:, iIdx)  = 2 .* chebA(:, iIdx-1)  + 2 .* x .* chebDA(:, iIdx-1)  - chebDA(:, iIdx-2);
                chebDDA(:, iIdx) = 4 .* chebDA(:, iIdx-1) + 2 .* x .* chebDDA(:, iIdx-1) - chebDDA(:, iIdx-2);
            end

            % assign to model properties
            model.A   = chebA;
            model.dA  = chebDA;
            model.ddA = chebDDA;

        end
        

        function state = updateC(model, state)

            state.c = model.A * state.coefC;
            
        end

        function state = updateW(model, state)

            state.w = model.A * state.coefW;

        end
        
        function state = updateDC(model, state)

            state.dC = model.dA * state.coefC;

        end
        
        function state = updateDW(model, state)

            state.dW = model.dA * state.coefW;

        end
        
        function state = updateDDC(model, state)

            state.ddC = model.ddA * state.coefC;

        end
        
        function state = updateDDW(model, state)

            state.ddW = model.ddA * state.coefW;
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



            state.eqC = massAccumC - dMobility .* dC .* dW - mobility .* ddW;

        end
        
        function state = updateEqW(model, state)
            % Residual form of : w + d2c/dx2 - F(c) = 0

            c    = state.c;
            ddC  = state.ddC;
            w    = state.w;
            
            F   = model.energyFunc(c);
            epsi = model.epsilon;

            state.eqW = w + epsi^2 .* ddC - F;
            
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

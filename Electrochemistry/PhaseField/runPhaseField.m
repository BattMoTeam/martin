%% read model parameters
%

filename = fullfile(battmoDir(), 'Electrochemistry', 'PhaseField', 'jsonfiles', 'phasefield.json');
jsonstruct = parseBattmoJson(filename);

% instantiate model (done also by setupPhaseFieldSimulation)
% inputparams = PhaseFieldInputParams(jsonstruct);
% model = PhaseField(inputparams);

%% Setup model, initial state
%

simsetup = setupPhaseFieldSimulation(jsonstruct);

%% Setup schedule
% in this case only the time steps are given, no source term
%

total = 100;
n     = 1000;                       
dt    = total / n;                  
dts   = rampupTimesteps(total, dt, 5);

clear flux
flux.functionFormat = 'tabulated';
flux.argumentList   = {'time'};
flux.dataX          = [0, 20, 21, 100];
flux.dataY          = [0.01, 0.01, -0.01, -0.01];

fluxfunc = setupFunction(flux);

control  = struct('src', @(t) fluxfunc(t));

step = struct('val', dts, 'control', ones(numel(dts), 1));
schedule = struct('control', control, 'step', step);

simsetup.schedule = schedule;

%% run simulation

simsetup.model.verbose = true;
% simsetup.run();



%% Visualization
[simResults, globvars, reports] = simsetup.run();
model = simsetup.model;
x     = model.chebyshevNodes(model.N);  
% Figure 1 : energy functions
c_test = linspace(0.01, 0.99, 200)';

% F'(c) used in eqW
f_prime = model.energyFunc(c_test);

% f_hom(c) double well shape
omega = model.omega;
kT    = model.kT;
f_hom = omega .* c_test .* (1 - c_test) + kT .* (c_test .* log(c_test) + (1 - c_test) .* log(1 - c_test));

figure;

subplot(1, 2, 1);
plot(c_test, f_prime, 'b-', 'LineWidth', 1.5);
xlabel('c');
ylabel("F'(c)");
title('Energy derivative (used in eqW)');
grid on;

subplot(1, 2, 2);
plot(c_test, f_hom, 'r-', 'LineWidth', 1.5);
xlabel('c');
ylabel('f_{hom}(c)');
title('Free energy : double well');
grid on;

% Figure 2 : concentration profiles
nPlots = numel(dts);
% nPlots = 100;
nPlots = min(nPlots, numel(simResults));  % in case fewer states are saved

figure;
hold on;

% initial state with black line
cInit = model.A * simsetup.initstate.coefC;
plot(x, cInit, 'k-', 'LineWidth', 2, 'DisplayName', 'initial');

% states 1 to nPlots with color gradient (blue = init, yellow = last)
cmap = parula(nPlots);
for iState = 1 : nPlots
    coefC = simResults{iState}.coefC;
    c     = model.A * coefC;
    plot(x, c, 'Color', cmap(iState, :), 'DisplayName', sprintf('state %d', iState));
end

if nPlots <= 10
    legend('show', 'Location', 'best');
end
cb = colorbar;
cb.Ticks      = [0, 1];
cb.TickLabels = {'Initial state', 'Final state'};
xlabel('x');
ylabel('c');
title(sprintf('Concentration profiles : initial + states 1 to %d', nPlots));
grid on;
hold off;

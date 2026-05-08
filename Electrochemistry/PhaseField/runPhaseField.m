%% read model parameters
%

filename = fullfile(battmoDir(), 'Electrochemistry', 'PhaseField', 'jsonfiles', 'phasefield.json');
jsonstruct = parseBattmoJson(filename);

%% Setup model, initial state
%

simsetup = setupPhaseFieldSimulation(jsonstruct);

%% Setup schedule
% in this case only the time steps are given, no source term
%

total = 1*hour;

n   = 100;
dt  = total/n;
dts = rampupTimesteps(total, dt, 5);

control  = struct('src', []);

step = struct('val', dts, 'control', ones(numel(dts), 1));
schedule = struct('control', control, 'step', step);

simsetup.schedule = schedule;

%% run simulation
%

simsetup.model.verbose = true;

simsetup.run();

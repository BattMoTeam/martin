filename = fullfile(battmoDir(), 'Electrochemistry', 'PhaseField', 'jsonfiles', 'phasefield.json');
jsonstruct = parseBattmoJson(filename);

simsetup = setupPhaseFieldSimulation(jsonstruct);

total = 1*hour;

n   = 100;
dt  = total/n;
dts = rampupTimesteps(total, dt, 5);

control  = struct('src', srcfunc);

step = struct('val', dts, 'control', ones(numel(dts), 1));
schedule = struct('control', control, 'step', step);

simsetup.schedule = schedule;

simsetup.model.verbose = true;

simsetup.run();

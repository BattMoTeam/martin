filename = fullfile(battmoDir(), 'Electrochemistry', 'PhaseField', 'jsonfiles', 'phasefield.json');
jsonstruct = parseBattmoJson(filename);

inputparams = PhaseFieldInputParams(jsonstruct);

model = PhaseField(inputparams);

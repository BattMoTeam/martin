%% material properties
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_lnmo_graphite.json');
jsonstruct_material = parseBattmoJson(jsonfilename);
jsonstruct_material.use_thermal = false;

filename = fullfile(battmoDir(), 'Electrochemistry', 'PhaseField', 'jsonfiles', 'phasefield.json');
jsonstruct_phasefield = parseBattmoJson(filename);

jsonstruct_material.PositiveElectrode.Coating.ActiveMaterial.diffusionModelType = 'phasefield';
jsonstruct_material.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion = jsonstruct_phasefield;

%% Control
% We load the json structure for the geometrical properties
jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
jsonstruct_control = parseBattmoJson(jsonfilename);

jsonstruct = mergeJsonStructs({jsonstruct_material, ...
                               jsonstruct_control});

model = setupModelFromJson(jsonstruct);

cgit = model.cgit;


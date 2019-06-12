% Quickly change the NKT parameters. Parameters to set are shown below

%% Parameters to set
varPWR = 100;       % Power level
varEMS = true;      % Emission state (true = on; false = off)
varLAM = 633;       % Center wavelenth on varia
varBWD = 3;         % Bandwidth on varia

%% Push params
fprintf('\n--- Setting parameters:\n')
fprintf('      Power      = %f\n', varPWR);
fprintf('      Emission   = %s\n', string(varEMS));
fprintf('      Wavelength = %f\n', varLAM);
fprintf('      Bandwidth  = %f\n', varBWD);

%-- Instantiate the NKT (will be "NKT" variable)
VFN_setUpNKT;   

%-- Set power
VFN_NKT_setPowerLevel(NKT, varPWR);

%-- Set emission state
VFN_NKT_setEmission(NKT, varEMS);

%-- Set wavelength on varia
VFN_NKT_setWvlRange(NKT, varLAM-varBWD/2, varLAM+varBWD/2);

%-- Close connection
VFN_cleanUpNKT;

fprintf('\n--- DONE\n');
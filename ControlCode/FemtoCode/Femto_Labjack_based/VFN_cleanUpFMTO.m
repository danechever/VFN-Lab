% Close the Labjack/Femto connection

% set gain to minimum (5)
FMTO.FMTO_scale = 5;
VFN_FMTO_setGain(FMTO);

%% Close the connection and populate the bench struct accordingly
py.labjack.ljm.close(FMTO.LJ);

%% the below commented lines are from HCST-R
%ljm.close doesn't return anything, and trying to read a value after
%closure will throw an error/halt the program, so just assume for now.
%bench.Femto.CONNECTED = false; 

%bench.Femto = rmfield(bench.Femto, 'LJ');

%disp('*** Labjack T7 disconnected. ***')
% Save backup bench object
%hcst_backUpBench(bench)
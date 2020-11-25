% Initialize the Labjack T7 and Femto for working with VFN.
% This function only opens the connection. (Modified from HCST-R library
% which also prints the device info in human-readable format, and prints 
% the current state of the Femto.)
%
% ALL Labjack functions will call the LJM Python library instead of using
% native Matlab functions, since there aren't any for Linux at this time.

%% Open the connection to the first T7 the computer can find and display device info
FMTO.LJ = py.labjack.ljm.openS("T7", "ANY", "ANY");

%% commented lines below are from HCST-R
%if(FMTO.s)
%    bench.Femto.CONNECTED = true;

%bench.Femto.LJ.info = py.labjack.ljm.getHandleInfo(bench.Femto.LJ.handle);

%sprintf("Opened a LabJack with Device type: %i, Connection type: %i,\nSerial number: %i, Port: %i,\nMax bytes per MB: %i", ...
%      bench.Femto.LJ.info{1}, bench.Femto.LJ.info{2}, bench.Femto.LJ.info{3}, bench.Femto.LJ.info{5}, bench.Femto.LJ.info{6})

%% Set DIO ports on the Labjack to be outputs
% py.labjack.ljm.eWriteName(bench.Femto.LJ.handle, "DIO_INHIBIT", int64(1));
  
%% Read the current Femto settings and display them
%hcst_readFemtoGain(bench);
%hcst_readFemtoOutput(bench);
%fprintf("Current Femto Voltage Signal:  %e V", bench.Femto.V)
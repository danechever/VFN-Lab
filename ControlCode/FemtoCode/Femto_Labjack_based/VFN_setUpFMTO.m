% Initialize the Labjack T7 and Femto for working with VFN.
% This function only opens the connection. (Modified from HCST-R library
% which also prints the device info in human-readable format, and prints 
% the current state of the Femto.)
%
% NOTE: user must provide the serial number of the device to connect to.
%       Remember to update this when a new device is connected!
%
% ALL Labjack functions will call the LJM Python library instead of using
% native Matlab functions, since there aren't any for Linux at this time.

%% USER-DEFINED Values
%-- Serial Number of device to connect to
FMTO.SSN = "470014583";

%-- AIN Port (which AIN# is used for the voltage/signal reads)
FMTO.AINPort = "AIN1";

%-- DIO Ports (which pins are used for each Femto Pin)
  % Order of Femto pins: {12, 11, 10, 14}
FMTO.DIOPorts = {"FIO2", "FIO1", "FIO0", "FIO3"};

%% Open the connection to the requested T7 SSN (using any connection type)
FMTO.LJ = py.labjack.ljm.openS("T7", "ANY", FMTO.SSN);

%% SET AIN1 port to +/-10V to prevent overvoltage
py.labjack.ljm.eWriteName(FMTO.LJ, strcat(FMTO.AINPort,"_RANGE"), int64(10));

%% Set DIO ports on the Labjack to be outputs
py.labjack.ljm.eWriteName(FMTO.LJ, "DIO_INHIBIT", int64(1));
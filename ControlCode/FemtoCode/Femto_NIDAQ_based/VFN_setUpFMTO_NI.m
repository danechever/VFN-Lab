%LUCI library
addpath('C:\Program Files (x86)\FEMTO\LUCI-10\Software\Matlab Libraries');

% Create Session
FMTO.s = daq.createSession('ni');

% Add analog input channel
addAnalogInputChannel(FMTO.s,'Dev1', 0, 'Voltage');

% Set Scan parameters
FMTO.s.Rate = FMTO.Nrate;             % [samples/second]
FMTO.s.NumberOfScans = FMTO.Nread;    % Number of samples at a given location

%call LUCI dll if not already called
if not(libisloaded('LUCI_10_x64'))
    loadlibrary('LUCI_10_x64', 'LUCI_10.h')
end
% Indexes all LUCI cables
calllib('LUCI_10_x64', 'EnumerateUsbDevices');
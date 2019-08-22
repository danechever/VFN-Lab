addpath(genpath('C:\Users\AOlab1\Documents\MATLAB\Add-Ons\Toolboxes\Zaber Device Control Toolbox\code'))

%LUCI library
addpath('C:\Program Files (x86)\FEMTO\LUCI-10\Software\Matlab Libraries');

% Create Session
s = daq.createSession('ni');

% Add analog input channel
addAnalogInputChannel(s,'Dev2', 0, 'Voltage');

% Set Scan parameters
s.Rate = Nrate;             % [samples/second]
s.NumberOfScans = Nread;    % Number of samples at a given location

%call LUCI dll if not already called
if not(libisloaded('LUCI_10_x64'))
    loadlibrary('LUCI_10_x64', 'LUCI_10.h')
end
% Indexes all LUCI cables
calllib('LUCI_10_x64', 'EnumerateUsbDevices');
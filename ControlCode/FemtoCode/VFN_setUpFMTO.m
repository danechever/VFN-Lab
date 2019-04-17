addpath(genpath('C:\Users\AOlab1\Documents\MATLAB\Add-Ons\Toolboxes\Zaber Device Control Toolbox\code'))

% Create Session
s = daq.createSession('ni');

% Add analog input channel
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');

% Set Scan parameters
s.Rate = Nrate;             % [samples/second]
s.NumberOfScans = Nread;    % Number of samples at a given location

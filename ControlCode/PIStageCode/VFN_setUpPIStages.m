% This will connect to the desired PI stages (via USB)
% It assumes only a single axis is connected per controller
% It assumes the single axis on all controllers is called '1'
%
% In the end, you will be left with two variables needed for control:
%   'Controller' = instance of controller object
%   'PIdevs'     = struct with instances of device (USB Connection) objects
%
% NOTE: This code no longer checks if the controllers are already connected
%   As such, it will try to reconnect to the controllers even if already
%   connected. This may cause issues.

%% Stage Serial Numbers
% Change these values if different controllers are connected
EXPECTED_AXES = [119063717; ...
                 119069544; ...
                 119063721]; 

%% Load PI MATLAB Driver GCS2
%  (if not already loaded)
addpath ( 'C:\Program Files (x86)\Physik Instrumente (PI)\Software Suite\MATLAB_Driver\' ); 

% The "controller" is actually the MATLAB driver. I kept the naming
% convention from PI.
if ( ~exist ( 'Controller', 'var' ) || ~isa ( Controller, 'PI_GCS_Controller' ) )
    Controller = PI_GCS_Controller ();
end


%% Find and validate available axes
% List all E-873 stages available via USB
devs = Controller.EnumerateUSB('E-873');

%-- Check that all found stages are expected
% Preallocate array with SNs
SNs = strings(numel(devs),1);
% Iterate through stages extracting the SN
for i = 1:numel(devs)
    % Find index where SN information starts
    ind = strfind(devs{i}, 'SN ');
    % Extract SN number and place into SNs array
    SNs(i) = num2str(sscanf(devs{i}(ind:end), 'SN %d'));
    
    % Check that this SN was expected
    if ~ismember(SNs(i), num2str(EXPECTED_AXES))
        error('An unrecgonized PI E-873 Controller was found: SN %s', SNs(i))
    end
end

% Delete unecessary variables
clear('ind', 'EXPECTED_AXES')


%% Start Connection
%-- USB Connection
% Iterate through discovered devices, connecting to each    
for i = 1:numel(SNs)
    % Create an entry in the PIdevs struct with the SN as the fieldname
    tag = sprintf('AX_%s', SNs(i));
    % Connect to controller 
    PIdevs.(tag) = Controller.ConnectUSB(char(SNs(i)));
end
% Get all the available axes names
axs = fieldnames(PIdevs);

clear('tag')

% NOTE: the following code was kept for reference but will require mods to work
%-- Other connection types
%  boolPIdeviceConnected = false; 
% if ( exist ( 'PIdevice', 'var' ) )
%     if ( PIdevice.IsConnected )
%         boolPIdeviceConnected = true; 
%     end
% end
% if ( ~(boolPIdeviceConnected ) )
%     % RS232
%     comPort = 1;          % Look at the Device Manager to get the rigth COM Port.
%     baudRate = 115200;    % Look at the manual to get the rigth bau rate for your controller.
%     PIdevice = Controller.ConnectRS232 ( comPort, baudRate );
%
%     % TCP/IP
%     ip = 'xxx.xxx.xx.xxx';  % Use "devicesTcpIp = Controller.EnumerateTCPIPDevices('')" to get all PI controller available on the network.
%     port = 50000;           % Is 50000 for almost all PI controllers
%     PIdevice = Controller.ConnectTCPIP ( ip, port ) ;
% end

% OPTIONAL: Query controller identification strings to confirm connection
for i = 1:numel(axs)
    fprintf('Connected to: %s',PIdevs.(axs{i}).qIDN());
end

%% Initialize PIdevs for use in MATLAB
for i = 1:numel(axs)
    PIdevs.(axs{i}) = PIdevs.(axs{i}).InitializeController();
end

%% Show connected stages
clear axs devs SNs i
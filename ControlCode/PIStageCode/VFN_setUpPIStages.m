% This will connect to the desired PI stages (via USB or Ethernet)
%  - User selects connection type at top of this script
% It assumes only a single axis is connected per controller
% It assumes the single axis on all controllers is called '1'
% It automatically imports MATLAB driver from the correct path for linux/pc
% 
% The USB mode automatically connects to all PI USB E-873 devices and
% checks that their SN matches a user-provided/expected SN
%
% The TCP mode connects only to user-provided IP addresses. It does not
% confirm any SN but still returns the PIdevs names with SN attached.
% NOTE: User must provide IPs in the requested text (see code below) 
%       This is to avoid putting our IPs on GIT
%
% In the end, you will be left with two variables needed for control:
%   'Controller' = instance of controller object
%   'PIdevs'     = struct with instances of device (USB Connection) objects
%
% NOTE: This code no longer checks if the controllers are already connected
%   As such, it will try to reconnect to the controllers even if already
%   connected. This may cause issues.

%% Connection Settings
% Change these values to change the connection type 

%-- Choose type ("usb" or "tcp")
contype = 'tcp';


%Make sure contype is lowercase
contype = lower(contype);

if strcmp(contype, 'usb')
    %-- USB
    % Provide serial numbers for devices since USB connection automatically
        % detects and connects to USB devices. 
    % Change these values if different controllers are connected
    USB_ExpectedSNs = [119063717; ...
                     119069544; ...
                     119063721];
elseif strcmp(contype, 'tcp')
    %-- TCP/IP
    % TO USER: enter the requested IP addresses into the file below.
    %          Also modify path to this file as needed
    flnm = '/home/vfndev/Documents/MATLAB/VFN-Lab/VFN_Config/VFN_PIStages_Q545_IPAdresses';
    
    % Read IP addresses from file and convert to char array
    IPs = table2array(readtable(flnm, 'NumHeaderLines',1));
else 
    %-- unrecognized contype
    error('Invalid connection type "contype"')
end

%% Load PI MATLAB Driver GCS2
%  (if not already loaded)

% Add path to the driver
if ispc
    % Windows Machine Detected. Use windows path:
    addpath ( 'C:\Program Files (x86)\Physik Instrumente (PI)\Software Suite\MATLAB_Driver\' ); 
elseif isunix
    % Linux Machine Detected. Use current path:
    pth = fileparts(mfilename('fullpath'));
    addpath(genpath([pth '/PI_Linux_Drivers']))
    %-- NOTE: this assumes the Matlab driver is in the directory which 
    %         holds this file. That is the case in the VFN-Lab git repo
end

clear pth

% The "controller" is actually the MATLAB driver. I kept the naming
% convention from PI.
if ( ~exist ( 'Controller', 'var' ) || ~isa ( Controller, 'PI_GCS_Controller' ) )
    Controller = PI_GCS_Controller ();
end


%% Find and validate available axes (in USB case)
if strcmp(contype, 'usb')
    %-- USB
    % List all E-873 stages available via USB
    devs = Controller.EnumerateUSB('E-873');
    %- Check that all found stages are expected
    % Preallocate array with SNs
    SNs = strings(numel(devs),1);
    % Iterate through stages extracting the SN
    for i = 1:numel(devs)
        % Find index where SN information starts
        ind = strfind(devs{i}, 'SN ');
        % Extract SN number and place into SNs array
        SNs(i) = num2str(sscanf(devs{i}(ind:end), 'SN %d'));

        % Check that this SN was expected
        if ~ismember(SNs(i), num2str(USB_ExpectedSNs))
            error('An unrecgonized PI E-873 Controller was found: SN %s', SNs(i))
        end
    end
end

% Delete unecessary variables
clear('ind', 'USB_ExpectedSNs')


%% Start Connection
if strcmp(contype, 'usb')
    %-- USB Connection
    % Iterate through discovered devices, connecting to each    
    for i = 1:numel(SNs)
        % Create an entry in the PIdevs struct with the SN as the fieldname
        tag = sprintf('AX_%s', SNs(i));
        % Connect to controller 
        PIdevs.(tag) = Controller.ConnectUSB(char(SNs(i)));
    end
elseif strcmp(contype, 'tcp')
    %-- TCP Connection
    % Preallocate array with SNs
    SNs = strings(numel(IPs),1);
    % Iterate through provided IP addresses, connecting to each
    for i = 1:numel(IPs)
        % Connect to controller (assuming port=50000 as typical for PI)
        tmpdev = Controller.ConnectTCPIP(IPs{i}, 50000);
        
        % Get serial number for PIdevs field name
        idn = strsplit(tmpdev.qIDN(), ',');
        SNs(i) = num2str(sscanf(idn{3}, '%d'));
        
        % Create an entry in the PIdevs struct with the SN as the fieldname
        tag = sprintf('AX_%s', SNs(i));
        % Connect to controller 
        PIdevs.(tag) = tmpdev;
    end
end

clear tag tmpdev idn IPs
    
% Get all the available axes names
axs = fieldnames(PIdevs);

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
clear axs devs SNs i contype
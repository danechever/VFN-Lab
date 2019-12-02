%% Comment:
% This will connect to the PI E-873 axis and initialize it.
% It assumes only a single axis is connected (this is the only option for
%       the E-873 anyway...)
% In the end, you will be left with three variables needed for control:
    % 1) 'axis'         = multidimensional cell containing the axis names for use by other scripts
    % 2) 'Controller'   = instance of controller object
    % 3) 'PIdevice'     = instance of device (USB connection) object

%% Load PI MATLAB Driver GCS2
%  (if not already loaded)
addpath ( 'C:\Users\Public\PI\PI_MATLAB_Driver_GCS2' ); % If you are still using XP, please look at the manual for the right path to include.

if ( ~exist ( 'Controller', 'var' ) || ~isa ( Controller, 'PI_GCS_Controller' ) )
    Controller = PI_GCS_Controller ();
end;


%% Start connection
%(if not already connected)

boolPIdeviceConnected = false; if ( exist ( 'PIdevice', 'var' ) ), if ( PIdevice.IsConnected ), boolPIdeviceConnected = true; end; end;

if ( ~(boolPIdeviceConnected ) )
% Please choose the connection type you need.    
%     % RS232
%     comPort = 1;          % Look at the Device Manager to get the rigth COM Port.
%     baudRate = 115200;    % Look at the manual to get the rigth bau rate for your controller.
%     PIdevice = Controller.ConnectRS232 ( comPort, baudRate );

%     % USB
    controllerSerialNumber = '0119017815';    % Use "devicesUsb = Controller.EnumerateUSB('')" to get all PI controller connected to you PC. Or look at the label of the case of your controller
    PIdevice = Controller.ConnectUSB ( controllerSerialNumber );

%     % TCP/IP
%     ip = 'xxx.xxx.xx.xxx';  % Use "devicesTcpIp = Controller.EnumerateTCPIPDevices('')" to get all PI controller available on the network.
%     port = 50000;           % Is 50000 for almost all PI controllers
%     PIdevice = Controller.ConnectTCPIP ( ip, port ) ;
end

% OPTIONAL: Query controller identification string
fprintf('\nConnected to: \n%s',PIdevice.qIDN());

% initialize PIdevice object for use in MATLAB
PIdevice = PIdevice.InitializeController ();


%% Show connected stages

% query names of all controller axes
availableAxes = PIdevice.qSAI_ALL;


% get Name of the stage connected to axis 1
% PIdevice.qCST ( '1' )
% Get name of stage connected to axis queried before
for iii = 1:length(availableAxes)
    fprintf('-- Stage [axis=type]: %s\n', PIdevice.qCST(char(availableAxes(iii))));
end

% UNEEDED: Show for all axes: which stage is connected to which axis
% for idx = 1 : 1%length ( availableAxes )
%     % qCST gets the name of the 
%     stageName = PIdevice.qCST ( availableAxes ( idx ) );
%     disp ( [ 'Axis ', availableAxes(idx), ': ', stageName ] );
% end

% if the stages listed are wrong use one of the following solutions:
% - PIMikroMove (recommended, but MS windows only)
% - PITerminal  (use VST, qCST, CST and WPA command as described in the Manual)
% - Use function "PI_ChangeConnectedStage"

% Define axis variable for single (first) axis in list
axis = availableAxes;

clear availableAxes boolPIdeviceConnected controllerSerialNumber iii
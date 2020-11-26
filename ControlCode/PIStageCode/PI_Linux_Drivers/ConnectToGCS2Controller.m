%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This sample demonstrates
% - Loading PI MATLAB Driver GCS2
% - Connecting to the controller
% - Retrieving controller Name
% - Closing the connection to the controller
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load PI MATLAB Driver GCS2
%  (if not already loaded)
addpath ( 'C:\Users\Public\PI\PI_MATLAB_Driver_GCS2' ); % If you are still using XP, please look at the manual for the right path to include.

if ( ~exist ( 'Controller', 'var' ) || ~isa ( Controller, 'PI_GCS_Controller' ) )
    Controller = PI_GCS_Controller ();
end;


%% Start connection
%(if not already connected)

% Test if PIdevice is connected
boolPIdeviceConnected = false; 
if ( exist ( 'PIdevice', 'var' ) )
    if ( PIdevice.IsConnected ), boolPIdeviceConnected = true; 
    end; 
end;


if ( ~(boolPIdeviceConnected ) )
% Please choose the connection type you need.    
    
%     % RS232
% On Linux use PI_ConnectRS232ByDevName.m instead
    comPort = 1;          % Look at the Device Manager to get the correct COM Port.
    baudRate = 115200;    % Look at the manual to get the correct bau rate for your controller.
    PIdevice = Controller.ConnectRS232 ( comPort, baudRate );


%     % USB
%     controllerSerialNumber = 'Enter valid serial number here, e.g. "123456789"';    % Use "devicesUsb = Controller.EnumerateUSB('')" to get all PI controller connected to you PC.
%                                                                                     % Or look at the label of the case of your controller
%     PIdevice = Controller.ConnectUSB ( controllerSerialNumber );


%    % TCP/IP
%     ip = 'xxx.xxx.xx.xxx';  % Use "devicesTcpIp = Controller.EnumerateTCPIPDevices('')" to get all PI controller available on the network.
%     port = 50000;           % Is 50000 for almost all PI controllers
%     PIdevice = Controller.ConnectTCPIP ( ip, port ) ;

end

% Query controller identification string
connectedControllerName = PIdevice.qIDN ()

% initialize PIdevice object for use in MATLAB
PIdevice = PIdevice.InitializeController ();


%% If you want to close the connection
PIdevice.CloseConnection;


%% if you want to unload the PI MATLAB Driver GCS2
Controller.Destroy;
clear Controller;
clear PIdevice;
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
    controllerSerialNumber = '0';    % Use "devicesUsb = Controller.EnumerateUSB('')" to get all PI controller connected to you PC.
%                                                                                     % Or look at the label of the case of your controller
    PIdevice = Controller.ConnectUSB ( controllerSerialNumber );


%     % TCP/IP
%     ip = 'xxx.xxx.xx.xxx';  % Use "devicesTcpIp = Controller.EnumerateTCPIPDevices('')" to get all PI controller available on the network.
%     port = 50000;           % Is 50000 for almost all PI controllers
%     PIdevice = Controller.ConnectTCPIP ( ip, port ) ;

end

% Query controller identification string
PIdevice.qIDN()

% initialize PIdevice object for use in MATLAB
PIdevice = PIdevice.InitializeController ();


%% Show connected stages

% query names of all controller axes
availableAxes = PIdevice.qSAI_ALL


% get Name of the stage connected to axis 1
PIdevice.qCST ( '1' )


% Show for all axes: which stage is connected to which axis
for idx = 1 : 1%length ( availableAxes )
    % qCST gets the name of the 
    stageName = PIdevice.qCST ( availableAxes ( idx ) );
    disp ( [ 'Axis ', availableAxes(idx), ': ', stageName ] );
end

% if the stages listed are wrong use one of the following solutions:
% - PIMikroMove (recommended, but MS windows only)
% - PITerminal  (use VST, qCST, CST and WPA command as described in the Manual)
% - Use function "PI_ChangeConnectedStage"


%% Startup Stage
axis = '1';

% switch servo on for axis
switchOn    = 1;
% switchOff   = 0;
PIdevice.SVO ( axis, switchOn );

% reference axis
PIdevice.FRF ( axis );  % find reference
bReferencing = 1;                           
disp ( 'Stage is referencing')
% wait for referencing to finish
while(0 ~= PIdevice.qFRF ( axis ) == 0 )                        
    pause(0.1);           
    [char ( 8 ) * ones( 1, 7 ), '.']
end       


%% Move Stage

% determine the allowed travel range of the stage
travelRangeMinimumPosition = PIdevice.qTMN ( axis )
travelRangeMaximumPosition = PIdevice.qTMX ( axis )

commandedPosition = rand ( 1 ) * ( travelRangeMaximumPosition - travelRangeMinimumPosition ) + travelRangeMinimumPosition
PIdevice.MOV ( axis, commandedPosition );

disp ( 'Stage is moving')
% wait for motion to stop
while(0 ~= PIdevice.IsMoving ( axis ) )
    pause ( 0.1 );
    [char ( 8 ) * ones( 1, 7 ), '.']
end

positonReached = PIdevice.qPOS(axis)


%% Call Sample Functions: PI_xxx ( PIdevice ... )
% This section shows some advanced features of your PIdevice. Please read
% the documentation first to get an idea of what the sample will show to
% avoid damage or unexpected behavior of your system.
% By the way working with separated functions (or scripts) is BEST PRACTICE
% to use the PI MATLAB driver and your PI device: Use this script to
% connect to your PI device and do the stage initialization (you may do
% some cleanup to get rid of unused code). Write your own script or
% function to implement your tasks - without doing the initialization stuff
% every time you run it.

% stageType = 'Enter valid stage type here, e.g. "M-123.XY42"';
% axisname = '1';
% PI_ChangeConnectedStage ( PIdevice, axis, stageType );

% PI_UseGcsCommandsDirectly ( PIdevice );

% PI_ListConnectedPiDevices ( Controller )

% PI_UseDataRecorder ( PIdevice, axis );


%% If you want to close the connection
PIdevice.CloseConnection ();

%% If you want to unload the dll and destroy the class object
Controller.Destroy ();
clear Controller;
clear PIdevice;
addpath(genpath('C:\Users\AOlab1\Documents\MATLAB\Add-Ons\Toolboxes\Zaber Device Control Toolbox\code'))

%% Axes Serial Numbers
% Change these values if the serial numbers on the axes are changed
AXIS_63 = 51463;    
AXIS_93 = 52893;    
AXIS_14 = 52714; 

%%  Prepare MATLAB serial port
% Instantiate the serial object
port = serial('COM3');

% Set default serial port properties for the binary protocol.
set(port, ...
    'BaudRate',     115200, ...
    'DataBits',     8, ...
    'FlowControl',  'none', ...
    'Parity',       'none', ...
    'StopBits',     1, ...
    'Terminator',   'CR/LF', ...
    'Timeout',      0.5);

% Disable timeout warning; Zaber lib intentionally times out for some reads
warning off MATLAB:serial:fgetl:unsuccessfulRead    % This is the warning for ASCII protocol
%warning off MATLAB:serial:fread:unsuccessfulRead   % This is the warning for Binary protocol

% Open the port.
fopen(port);

%% Instantiate the Zaber Classes
try
    % Create AsciiProtocol object (only 1 since axes are daisy-chained)
    protocol = Zaber.AsciiProtocol(port);
    
    devs = protocol.finddevices();
    
    % Use device serial numbers to match axes
    for i = 1:length(devs)
        serNum = devs(i).request('get', 'system.serial').Data;
        if serNum == AXIS_63
            ax63 = devs(i);
        elseif serNum == AXIS_93
            ax93 = devs(i);
        elseif serNum == AXIS_14
            ax14 = devs(i);
        else
            error('An unrecognized zaber was found\n  Serial Numer: %d',...
                serNum)
        end
    end
    
%% Home the axes if they have not been homed
    if exist('ax14', 'var')
        if ~ax14.request('get', 'limit.home.triggered').Data
            ax14.home;
            fprintf('ax14 Homed\n')
        end
    else
        fprintf('ax14 Not Present\n')
    end
        
    if exist('ax63', 'var')
        if ~ax63.request('get', 'limit.home.triggered').Data
            ax63.home;
            fprintf('ax63 Homed\n')
        end
    else
        fprintf('ax63 Not Present\n')
    end
    
     if exist('ax93', 'var')
        if ~ax93.request('get', 'limit.home.triggered').Data
            ax93.home;
            fprintf('ax93 Homed\n')
        end
    else
        fprintf('ax93 Not Present\n')
    end   
    
catch exception
    % Close port if an error occurs, otherwise it remains locked
    fclose(port);
    rethrow(exception);
end

clear('AXIS*', 'protocol', 'serNum', 'devs', 'i')

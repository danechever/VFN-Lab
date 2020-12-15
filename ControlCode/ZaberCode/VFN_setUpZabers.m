% This will connect to the desired Zaber stages
% It assumes all Zabers are daisy-chained and connected to a single USB
%
% NOTE: This script will home any Zabers that are unreferenced
%
% Note: the "port" variable is now contained in the "Zabs" struct
% 
% In the end, you will be left with a single variable needed for control:
%   'Zabs'  = struct with instances of device objects
%

%% Imports and setup
% Import Zaber objects
import zaber.motion.Library;
import zaber.motion.ascii.Connection;
import zaber.motion.ascii.AxisSettings


%-- Allow library to download the latest Database 
% Stores device database so that it doesn't need to be pulled from internet
% again in the future. More details here: 
% - https://www.zaber.com/software/docs/motion-library/ascii/howtos/device_db/
Library.enableDeviceDbStore();

%% Axes Serial Numbers
% Change these values if the serial numbers on the axes are changed
    % This feature is technically not needed but I decided to add it since
    % it forces the user to be conscious of which zaber devices are
    % connected. Like this they are aware of any numbering issues and such.
EXPECTED_AXES = [51463, ...
                 52893, ...
                 52714, ...
                 59214, ...
                 60707]; 

%%  Open serial port/connection
% Use the zaber library to open the port
Zabs.port = Connection.openSerialPort('/dev/ttyZabers');

%% Instantiate the Zaber Classes
try
    % Find all connected devices
    devs = Zabs.port.detectDevices();
    
    % Use device serial numbers to match axes and create instance in struct
    for i = 1:devs.length
        serNum = devs(i).getSerialNumber;
        if ismember(serNum, EXPECTED_AXES)
            % Device was expected; add it to the output struct
            tag = sprintf('AX_%d', serNum);
            Zabs.(tag) = devs(i).getAxis(1);
        else
            error('An unrecognized zaber was found\n Serial Number :%d',...
                serNum)
        end
    end
    
    %% Home each axis if it has not been homed
    % Get all elements in the struct so that we can iterate through them
    axs = fieldnames(Zabs);
    
    % Iterate through axes
    for i = 1:numel(axs)
        if ~isa(Zabs.(axs{i}), 'zaber.motion.ascii.Axis')
            % Skip non-zaber object fields (ie. port object)
            continue
        end
        if ~Zabs.(axs{i}).getSettings().get('limit.home.triggered')
            % Axis was not homed so home it
            Zabs.(axs{i}).home
            fprintf('%s Homed\n', axs{i})
        end
    end
    
catch exception
    % Close port if an error occurs, otherwise it remains locked
    Zabs.port.close();
    fprintf('An error occurred, port was closed\n');
    rethrow(exception);
end

clear('EXPECTED_AXES', 'serNum', 'devs', 'axs', 'tag', 'i')

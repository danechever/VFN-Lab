%Set up the Thorlabs PM100D for use
%- Connects to power meter
%- Creates a "PM" variable containing pertinent object information
%- Attempts to automatically connect to first PM100D found by system
%- Commented-out section allows manual entry of PM100D info
%
%* Attempts ot connect automatically on Windows or Linux
%
%This uses the ThorlabsPM100 python library
%- To install: https://github.com/clade/ThorlabsPM100
%- Relies on PyVisa: https://pyvisa.readthedocs.io/en/stable/
%- Also NI-VISA: https://www.ni.com/en-us/support/downloads/drivers/download.ni-visa.html#329456
%
% See also VFN_PM100D_read, VFN_cleanUpPM100D
%
%  author: D. Echeverri (dechever@caltech.edu)
%  created: 02/25/2020

if ispc
    %-- Windows PC detected. Act accordingly
    % Define manufacturer ID and product ID for auto-device detect
    PM_IDS = '0x1313::0x8078';

    % Instantiate PyVisa object
    PM.rm = py.visa.ResourceManager();

    %-- Connect to power meter (find first PM100D device automatically)
    % List all devices visible to PyVisa
    PM.devs = PM.rm.list_resources();
    % Iterate through list looking for first PM100D device
    for devs_i = 1:py.len(PM.devs)
        if strfind(char(PM.devs{devs_i}), PM_IDS)
            % Device with matching Manufacturer and Product ID was found
            PM.dev = char(PM.devs{devs_i});
            break
        end
        if devs_i == py.len(PM.devs)
            % Looked at all devices; none matched our ID string
            error('No PM100D was auto-detected. Try manually connecting (see setUp code)')
        end
    end
    % Instantiate instrument object and connect to device
    PM.inst = PM.rm.open_resource(PM.dev, pyargs('timeout',1));

    % %-- Manual way to connect to PM100D
    % %- Set the Serial Number (P0022602 in this exmple) to match your device
    % try
    %     PM.inst = PM.rm.open_resource('USB0::0x1313::0x8078::P0022602::0::INSTR',pyargs('timeout',1));
    % catch pyERR
    %     disp(getReport(pyERR,'basic'))
    %     error('Error connecting to device. Check serial number and that device is connected')
    % end
elseif isunix
    %-- Linux Machine detected. Act accordingly
    % Connect using USBTMC to default device (/dev/usbtmc0)
    PM.inst = py.ThorlabsPM100.USBTMC();
    
	% %-- If default USBTMC path fails, provide path here
	% PM.inst = py.ThorlabsPM100.USBTMC(pyargs('device','/dev/usbtmc0'));
    
    % Alias PM.rm to PM.inst so that cleanUp func. will work correctly
    PM.rm = PM.inst;
end

% Instantiate device main comms object
PM.PM = py.ThorlabsPM100.ThorlabsPM100(PM.inst);

% Delete unecessary values
clear('devs_i', 'PM_IDS')
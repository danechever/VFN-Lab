%Set up TWO Thorlabs PM100D power meters for simultaneous use
%- Connects to both power meters
%- Creates "PM" variables containing pertinent object information
%
%* ONLY TESTED ON LINUX
%
%This uses the ThorlabsPM100 python library
%- To install: https://github.com/clade/ThorlabsPM100
%
% See also VFN_PM100D_read, VFN_cleanUpPM100D
%
%  author: D. Echeverri (dechever@caltech.edu)
%  created: 07/12/2021

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
    % PM.inst = py.ThorlabsPM100.USBTMC();
    
	%-- If default USBTMC path fails, provide path here
    % Provide the SNs for each of the two redPMs to use
    normPM_SN = 'P0022602';
    nullPM_SN = 'P0015519';
    % Provide system path to both redPMs
    devPth1 = '/dev/ttyPM100D602';
    devPth2 = '/dev/ttyPM100D519';
    
    % Connect to the first PM
	instTMP = py.ThorlabsPM100.USBTMC(pyargs('device',devPth1));
    if count(string(instTMP.query('*IDN?')),normPM_SN)
        % The PM on devPth1 is the norm PM, save accordingly
        normPM.inst = instTMP;
    elseif count(string(instTMP.query('*IDN?')),nullPM_SN)
        % The PM on devPth1 is the nulling PM, save accordingly
        nullPM.inst = instTMP;
    else
        error('The redPM serial number could not be matched to a user-provided value')
    end
    
    % Connect to the second PM
	instTMP = py.ThorlabsPM100.USBTMC(pyargs('device',devPth2));
    if count(string(instTMP.query('*IDN?')),normPM_SN)
        % The PM on devPth2 is the norm PM, save accordingly
        normPM.inst = instTMP;
    elseif count(string(instTMP.query('*IDN?')),nullPM_SN)
        % The PM on devPth2 is the nulling PM, save accordingly
        nullPM.inst = instTMP;
    else
        error('The redPM serial number could not be matched to a user-provided value')
    end
end

% Instantiate device main comms objects
normPM.PM = py.ThorlabsPM100.ThorlabsPM100(normPM.inst);
nullPM.PM = py.ThorlabsPM100.ThorlabsPM100(nullPM.inst);

% Delete unecessary values
clear('devs_i', 'PM_IDS','normPM_SN','nullPM_SN','devPth1','devPth2','instTMP')
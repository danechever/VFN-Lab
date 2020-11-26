function devices = EnumerateTCPIPDevicesAsArray(c, szFilter)
% Get description strings of all PI devices on the network.
%
%   SYNTAX 
%       devices = ENUMERATETCPIPDEVICES()
%       devices = ENUMERATETCPIPDEVICES(szFilter)
%
%   INPUT 
%       szFilter           Filter on a special char sequence, e.g. 'E-042'
%       (char)             to list only the E-042 controllers available. 
%                           Leave empty to list all pi controllers available. 
%                           
%   OUTPUT 
%       devices            List of PI devices.
%       (cell array)
%
%   EXAMPLE
%       devices = ENUMERATETCPIPDEVICES()
%       devices = ENUMERATETCPIPDEVICES('E-042')
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

warning('"EnumerateTCPIPDevicesAsArray()" is deprecated und may not be supported in future versions of PI MATLAB Driver. Use "EnumerateTCPIPDevices()" instead.');

%% Settings


%% Function arguments handle
if nargin<2; szFilter=''; end

%% Check input parameters


%%  Program start

functionName = [ c.libalias, '_', mfilename];

if(any(strcmp(functionName,c.dllfunctions)))
    devices = blanks(10001);
    try
        [bRet,devices] = calllib(c.libalias,functionName,devices,10000,szFilter);
        if (bRet==0)
            devices = '';% no devices found -> return empty string.
        end
        
        devices = regexp(devices, '\n', 'split');
        devices = devices';
        
        % Usually the last entry is empty
        if (isempty(devices{end}))
            devices(end,:) = [];
        end
        
     catch
        rethrow(lasterror);
    end
else
    error('%s not found',functionName);
end
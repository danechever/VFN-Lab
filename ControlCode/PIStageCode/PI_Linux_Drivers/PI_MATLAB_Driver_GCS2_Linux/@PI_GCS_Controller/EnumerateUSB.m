function devices = EnumerateUSB(c, szFilter)
% Get description strings of all PI devices on the network.
%
%   SYNTAX 
%       devices = ENUMERATEUSB()
%       devices = ENUMERATEUSB(szFilter)
%
%   DESCRIPTION 
%   This function lists all PI devices connected to the PC. Make shure
%   you index the returned cell array "devices" with curly brackets to get
%   the string: stringOfFirstController = devices{1}
%    
%   INPUT 
%       szFilter           Filter on a special char sequence, e.g. 'E-042'
%       (char)             to list only the E-042 controllers available. 
%                          Leave empty to list all pi controllers available. 
%                           
%   OUTPUT 
%       devices            List of PI devices.
%       (cell array)
%
%   EXAMPLE
%       devices = ENUMERATEUSB()
%       devices = ENUMERATEUSB('E-042')
% 
%
%   COPYRIGHT © PHYSIKINSTRUMENTE (PI) GMBH U. CO. KG, support-software@pi.ws


%% Settings


%% Function arguments handle
if nargin<2; szFilter=''; end

%% Check input parameters


%%  Program start

functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
    szDevices = blanks(10001);
    try
        [bRet,szDevices] = calllib(c.libalias,functionName,szDevices,10000,szFilter);
        if (bRet==0)
            szDevices = '';% no devices found -> return empty string.
        end

        devices = regexp(szDevices, '\n', 'split');
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

function [c, iNumDevices]  = OpenUSBDaisyChain(c,serialnumber)

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

% Create variables for C interface
iNumDevices = 0;
piNumDevices = libpointer('int32Ptr',iNumDevices);
szDCDevices = blanks(2001);
try
    % call C interface
    % second returning argument would be "serialnumber" which is not needed. Hence is "~"
    [ c.DC_ID, ~, iNumDevices, szDCDevices]  = calllib(c.libalias,'PI_OpenUSBDaisyChain',serialnumber, piNumDevices,szDCDevices,2000);
    
    if (c.DC_ID<0)
        iError = GetError(c);
        szDesc = TranslateError(c,iError);
        error(szDesc);
    end
    
catch
    rethrow(lasterror);
end

% Create list of available Daisy Chain Devices
daisyChainDevicesAsCellArray = regexp(szDCDevices, '\n', 'split');

for address = 1 : size ( daisyChainDevicesAsCellArray , 2)
    if ( ~isempty ( daisyChainDevicesAsCellArray {address} ) )
        c.ConnectedDaisyChainDevices{address} = [sprintf('%02d: ', address), daisyChainDevicesAsCellArray{address}];
    end
end
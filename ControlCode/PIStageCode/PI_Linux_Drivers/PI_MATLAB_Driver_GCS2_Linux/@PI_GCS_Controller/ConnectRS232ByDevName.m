function c = ConnectRS232ByDevName(c, deviceName, baudRate)
% Open an RS-232 interface to a controller with Linux.
% 
%   SYNTAX 
%       E042  = PI_ConnectRS232ByDevName(deviceName, baudRate)
% 
%   DESCRIPTION 
%   The call also sets the baud rate on the controller side.
%    
%   INPUT 
%       deviceName              device interface name for RS232 connection
%       (char array)            (Linux).
%   
%       baudRate                to use.
%       (double)                
%                           
%   OUTPUT 
%       E042                    New controller object to communicate with a
%       (PI_GCS_Controller)     specific controller.
%       
%   EXAMPLE
%       E042  = PI_ConnectRS232ByDevName(deviceName, baudRate)
%
%   E042 is an alias for an arbitrary name of a controller object.

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% only load dll if it wasn't loaded before
if(~libisloaded(c.libalias))
	c = LoadGCSDLL (c);
end
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	try
		[c.ID] = calllib(c.libalias, functionName, deviceName, baudRate);
		if(c.ID<0)
			iError = GetError(c);
			szDesc = TranslateError (c, iError);
			error(szDesc);
		end
    catch ME
        error (ME.message);
	end
else
	error('%s not found',functionName);
end

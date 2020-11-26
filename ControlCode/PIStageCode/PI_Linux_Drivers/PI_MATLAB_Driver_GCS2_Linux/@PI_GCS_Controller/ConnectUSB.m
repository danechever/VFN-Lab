function [c] = ConnectUSB(c,szIdentifier)
%function [c] = ConnectUSB(c,szIdentifier)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	if(nargin<2),    szIdentifier = blanks(5000); end,
	try
		[c.ID, ~] = calllib(c.libalias, functionName, szIdentifier);
		if (c.ID == -1)
			error ( 'Interface could not be opened or no controller is responding. Check if controller is switched on and properly connected.' );
		end
	catch
		rethrow ( lasterror );
	end
else
	error('%s not found',functionName);
end

InitializeController(c);
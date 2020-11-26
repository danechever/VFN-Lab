function c = CloseDaisyChain(c)
%function CloseDaisyChain(c)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.DC_ID<0), error('The daisy chain is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	try
		calllib(c.libalias,functionName,c.DC_ID);
		c = SetDefaults(c);
	catch
		rethrow(lasterror);
	end
else
	error('%s not found',functionName);
end
c.DC_ID = -1;
c.ConnectedDaisyChainDevices = '';

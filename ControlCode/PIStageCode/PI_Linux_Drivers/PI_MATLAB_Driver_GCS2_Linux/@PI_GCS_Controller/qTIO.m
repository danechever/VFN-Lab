function [iNInputs,iNOutputs] = qTIO(c)
%function [iNInputs,iNOutputs] = qTIO(c)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	iNInputs = zeros(1);
	iNOutputs = zeros(1);
	piNInputs = libpointer('int32Ptr',iNInputs);
	piNOutputs = libpointer('int32Ptr',iNOutputs);
	try
		[bRet,iNInputs,iNOutputs] = calllib(c.libalias,functionName,c.ID,piNInputs,piNOutputs);
		if(bRet==0)
			iError = GetError(c);
			szDesc = TranslateError(c,iError);
			error(szDesc);
		end
	catch
		rethrow(lasterror);
	end
else
	error('%s not found',functionName);
end

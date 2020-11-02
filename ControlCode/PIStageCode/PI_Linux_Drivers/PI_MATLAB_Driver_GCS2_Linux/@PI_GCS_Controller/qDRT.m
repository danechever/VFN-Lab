function [iTriggerOptions] = qDRT(c,iChannels)
%function [iTriggerOptions] = qDRT(c,iChannels)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piChannels = libpointer('int32Ptr',iChannels);
	iTriggerOptions = zeros(size(iChannels));
	piTriggerOptions = libpointer('int32Ptr',iTriggerOptions);
	nValues = length(iTriggerOptions);
	szBuffer = blanks(100);
	try
		[bRet,iChannels,iTriggerOptions,szBuffer] = calllib(c.libalias,functionName,c.ID,piChannels,piTriggerOptions,szBuffer,nValues,99);
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

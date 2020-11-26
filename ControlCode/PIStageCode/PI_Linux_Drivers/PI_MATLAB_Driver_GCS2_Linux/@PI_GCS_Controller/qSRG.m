function [iOutValues1] = qSRG(c,szAxes,iInValues1)
%function [iOutValues1] = qSRG(c,szAxes,iInValues1)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piInValues1 = libpointer('int32Ptr',iInValues1);
	iOutValues1 = zeros(size(iInValues1));
	piOutValues1 = libpointer('int32Ptr',iOutValues1);
	try
		[bRet,szAxes,iInValues1,iOutValues1] = calllib(c.libalias,functionName,c.ID,szAxes,piInValues1,piOutValues1);
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

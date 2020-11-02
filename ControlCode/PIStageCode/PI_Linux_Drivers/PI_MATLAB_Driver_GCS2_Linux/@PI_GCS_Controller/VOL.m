function [bRet] = VOL(c,iInValues1,dInValues1)
%function [bRet] = VOL(c,iInValues1,dInValues1)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piInValues1 = libpointer('int32Ptr',iInValues1);
	pdInValues1 = libpointer('doublePtr',dInValues1);
	nValues = length(iValues1);
	try
		[bRet,iInValues1,dInValues1] = calllib(c.libalias,functionName,c.ID,piInValues1,pdInValues1,nValues);
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

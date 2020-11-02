function [bRet] = SPA(c,szAxes,uiInParamIDs,dInValues,szInString)
%function [bRet] = SPA(c,szAxes,uiInParamIDs,dInValues,szInString)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	if(nargin<5),szInString = '';end,
	puiInParIDs = libpointer('uint32Ptr',uiInParamIDs);
	pdInValues = libpointer('doublePtr',dInValues);
	try
		[bRet,szAxes,uiInParamIDs,dInValues,szInString] = calllib(c.libalias,functionName,c.ID,szAxes,puiInParIDs,pdInValues,szInString);
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

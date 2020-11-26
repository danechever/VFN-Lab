function [dOutValues1,szAnswer] = qSPA(c,szAxes,uiInParamIDs)
%function [dOutValues1,szAnswer] = qSPA(c,szAxes,uiInParamIDs)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	puiInPars = libpointer('uint32Ptr',uiInParamIDs);
	dOutValues1 = zeros(size(uiInParamIDs));
	pdOutValues1 = libpointer('doublePtr',dOutValues1);
	szAnswer = blanks(1001);
	try
		[bRet,szAxes,uiInParamIDs,dOutValues1,szAnswer] = calllib(c.libalias,functionName,c.ID,szAxes,puiInPars,pdOutValues1,szAnswer,1000);
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

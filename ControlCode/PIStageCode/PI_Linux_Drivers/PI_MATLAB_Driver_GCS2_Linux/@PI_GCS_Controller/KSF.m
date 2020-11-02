function bRet = KSF(c, szCoordinateSystemName)
% function bRet = KSF(c, szCoordinateSystemName)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

functionName = [ c.libalias, '_', mfilename];

if(c.ID < 0)
    error('The controller is not connected')
end;

if(~any( strcmp(functionName, c.dllfunctions )))
    error('%s not found',functionName);
else
	try
		[bRet] = calllib(c.libalias, functionName, c.ID, szCoordinateSystemName);
		if(bRet==0)
			iError = GetError(c);
			szDesc = TranslateError(c,iError);
			error(szDesc);
		end
	catch
		rethrow(lasterror);
    end
end
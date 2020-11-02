function szAnswer = qHDR(c)
% function szAnswer = qHDR(c)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	szAnswer = blanks(8001);
	try
		[bRet,szAnswer] = calllib(c.libalias,functionName,c.ID,szAnswer,8000);
		if(bRet==0)
			iError = GetError(c);
			if(iError == -5)
				szAnswer = blanks(18001);
				[bRet,szAnswer] = calllib(c.libalias,functionName,c.ID,szAnswer,18000);
				if(bRet==0)
					iError = GetError(c);
					szDesc = TranslateError(c,iError);
					error(szDesc);
				end
			else
				szDesc = TranslateError(c,iError);
				error(szDesc);
			end
		end
	catch
		rethrow(lasterror);
	end
else
	error('%s not found',functionName);
end

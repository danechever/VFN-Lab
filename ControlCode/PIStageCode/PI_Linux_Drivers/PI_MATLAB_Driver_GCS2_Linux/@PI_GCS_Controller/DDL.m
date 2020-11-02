function bRet = DDL(c,TableNr,StartIndex,Values)
%function bRet = DDL(c,iInValues1,iInValues2,dInValues1)

% This code is provided by Physik Instrumente(PI) GmbH&Co.KG
% You may alter it corresponding to your needs
% Comments and Corrections are very welcome
% Please contact us by mailing to support-software@pi.ws

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	pdValues = libpointer('doublePtr',Values);
	nValues = length(Values);
	try
		[bRet] = calllib(c.libalias,functionName,c.ID,TableNr,StartIndex,nValues,pdValues);
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

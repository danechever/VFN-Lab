function [szSourceIds,iRecordOptions] = qDRC(c,iRecordTables)
%function [szSourceIds,iRecordOptions] = qDRC(c,iRecordTables)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piRecordTables = libpointer('int32Ptr',iRecordTables);
	iRecordOptions = zeros(size(iRecordTables));
	piRecordOptions = libpointer('int32Ptr',iRecordOptions);
	nValues = length(iRecordTables);
	szSourceIds = blanks(1001);
	try
		[bRet,iRecordTables,szSourceIds,iRecordOptions] = calllib(c.libalias,functionName,c.ID,piRecordTables,szSourceIds,piRecordOptions,1000,nValues);
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

function [dOutValues1,szHeader] = qJLT(c,iJoystickIDs,iJoystickAxisIDs,iOffset,iNumberOfValues)
%function [dOutValues1,szHeader] = qJLT(c,iJoystickIDs,iJoystickAxisIDs,iOffset,iNumberOfValues)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piJoystickIDs = libpointer('int32Ptr',iJoystickIDs);
	piJoystickAxisIDs = libpointer('int32Ptr',iJoystickAxisIDs);
	szHeader = blanks(1000);
	ppdData = libpointer('doublePtr');
	nTables = length(iJoystickIDs);
	try
		[bRet,iJoystickIDs,iJoystickAxisIDs,ppdData,szHeader] = calllib(c.libalias,functionName,c.ID,piJoystickIDs,piJoystickAxisIDs,nTables,iOffset,iNumberOfValues,ppdData,szHeader,999);
		if(bRet==0)
			iError = GetError(c);
			szDesc = TranslateError(c,iError);
			error(szDesc);return;
		end
	headerlines = regexp(szHeader,'/n','split');
	for n = 1:length(headerlines)
		if(~isempty(strfind(headerlines{n},'DIM'))),nTables = str2num(headerlines{n}(strfind(headerlines{n},'=')+1:end));end;
		if(~isempty(strfind(headerlines{n},'NDATA'))),iNumberOfValues = str2num(headerlines{n}(strfind(headerlines{n},'=')+1:end));end;
	end
		i = 0;
		while(i<(nTables*iNumberOfValues))
			pause(0.1);
			i =  GetAsyncBufferIndex(c);
		end
		setdatatype(ppdData,'doublePtr',nTables,iNumberOfValues);
		dOutValues1 = ppdData.Value';
	catch
		rethrow(lasterror);
	end
else
	error('%s not found',functionName);
end

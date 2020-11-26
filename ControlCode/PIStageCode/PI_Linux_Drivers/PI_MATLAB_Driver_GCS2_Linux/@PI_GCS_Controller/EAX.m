function bRet = EAX(c,szAxes,iValues)


if(c.ID<0), error('The controller is not connected'),end;
AssertValueVectorSize( c, szAxes, iValues )
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piValues = libpointer('int32Ptr',iValues);
	try
		[bRet,szAxes,iValues] = calllib(c.libalias,functionName,c.ID,szAxes,piValues);
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

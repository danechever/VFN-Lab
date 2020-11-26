function bRet = WAV_SIN_P(c,WaveTableId, OffsetOfFirstPointInWaveTable, NumberOfPoints, AddAppendWave, CenterPointOfWave, AmplitudeOfWave, OffsetOfWave, SegmentLength)
% function bRet = WAV_SIN_P(c,WaveTableId, OffsetOfFirstPointInWaveTable, NumberOfPoints, AddAppendWave, CenterPointOfWave, AmplitudeOfWave, OffsetOfWave, SegmentLength)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName, c.dllfunctions)))
% 	pWaveTableId = libpointer('int32Ptr', WaveTableId);
%     pOffsetOfFirstPointInWaveTable = libpointer('int32Ptr', OffsetOfFirstPointInWaveTable); 
%     pNumberOfPoints = libpointer('int32Ptr', NumberOfPoints);
%     pAddAppendWave = libpointer('int32Ptr', AddAppendWave);
%     pCenterPointOfWave = libpointer('int32Ptr', CenterPointOfWave); 
%     pAmplitudeOfWave = libpointer('doublePtr', AmplitudeOfWave);
%     pOffsetOfWave = libpointer('doublePtr', OffsetOfWave);
%     pSegmentLength = libpointer('int32Ptr', SegmentLength);
    try
		[bRet] = calllib(c.libalias, functionName, c.ID, WaveTableId, OffsetOfFirstPointInWaveTable, NumberOfPoints, AddAppendWave, CenterPointOfWave, AmplitudeOfWave, OffsetOfWave, SegmentLength);
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

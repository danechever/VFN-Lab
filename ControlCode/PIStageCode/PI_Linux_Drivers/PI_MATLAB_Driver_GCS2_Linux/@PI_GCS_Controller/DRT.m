function DRT (c, recordTableId, triggerSource, triggerSourceValue)
% Configure Data Recorder Trigger.
% 
%   SYNTAX 
%   DRT (recordTableId, triggerSource, triggerSourceValue)
% 
%   DESCRIPTION 
%   Derfines a trigger source for data recorder. If the condition of the
%   trigger is fulfiled, the data recorder will start recording.
%    
%   INPUT 
%       recordTableId           Usualy "0". For further options please read
%       (integer)               the controller manual.
%   
%       triggerSource           Defines the condition which triggers the
%       (integer)               data recorder. Read the controller manual
%                               or use function "qHDR()" to get all
%                               available options for triggerSource.
%                               Frequently used options are (may not work
%                               for all controllers): 
%                               0 = default setting 
%                               1 = any command changing position (e.g. MOV) 
%                               2 = next command 
%     
%       triggerSourceValue      Additional parameter for special
%       (string)                configuration. Set to "'0'" if not needed,
%                               otherwise read the manual.
%       
%   EXAMPLE
%       This example may not work for all controllers. Usually it means:
%       Start recording if any command which changes target position is set (-> 1).
%       PI_Controller.DRT (0, 1, '0')

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	piValues1 = libpointer('int32Ptr',recordTableId);
	piValues2 = libpointer('int32Ptr',triggerSource);
	iValues3 = length(triggerSource);
	try
		[bRet,piValues1,piValues2,triggerSourceValue] = calllib(c.libalias,functionName,c.ID,piValues1,piValues2,triggerSourceValue,iValues3);
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

function DRC (c, recordTableId, recordSourceId, recordOptionId)
% Configure Data Recorder.
% 
%   SYNTAX 
%   DRC (recordTableId, recordSourceId, recordOptionId)
% 
%   DESCRIPTION 
%   This function configures the controller's data recorder. Usually a
%   specific kind of data (e.g. commanded or actual position) from a
%   specific record source (e.g. axis or channel) is recorded. The will be
%   stored in the table specified.
%    
%   INPUT 
%       recordTableId           ID of the table the data will be stored in.
%       (integer)               Needed for function "qDRR()".
%   
%       recordSourceId          Can be axis or channel name.
%       (string)                
%     
%       recordOptionId          Describes what kind of data will be
%       (integer)               recorded. A full list is provided in the
%                               controller manual or by function "qHDR()".
%                               Frequently used options are:
%                               1 = Commanded position of axis
%                               2 = Actual position of axis
%       
%   EXAMPLE
%       Record current position (-> 2) of axis/channel 3 (-> '3') into record table 1 (-> 1).
%       PI_Controller.DRC (1, '3', 2)

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	pRecordTableId  = libpointer('int32Ptr',recordTableId);
	pRecordOptionId = libpointer('int32Ptr',recordOptionId);
	try
		[bRet,pRecordTableId,recordSourceId,pRecordOptionId] = calllib(c.libalias,functionName,c.ID,pRecordTableId,recordSourceId,pRecordOptionId);
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
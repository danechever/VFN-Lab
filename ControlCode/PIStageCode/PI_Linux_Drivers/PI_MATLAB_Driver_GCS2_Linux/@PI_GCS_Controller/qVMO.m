function movePossible  = qVMO(c, axes, positions)
% Virtual move.
% 
%   SYNTAX 
%       movePossible  = E042.qVMO(axes, positons)
% 
%   DESCRIPTION 
%   Checks whether the moving platform of the Hexapod and additional single
%   axes can approach a specified position from the current position. 
%   Returns a boolean for all axes querried. qVMO(...) does not trigger any
%   motion. 
%    
%   INPUT 
%       axes                    String with axes. Write {'a', 'b'} to
%       (cell array of strings) use axes "a" and "b".
%   
%       positions               Target positions.
%       (double)                
%                           
%   OUTPUT 
%       movePossible            Indicates whether the 
%       (logical)               moving platform can approach the position 
%                               resulting from the given target position values:
%                               false = specified position cannot be approached
%                               true = specified position can be approached
%       
%   EXAMPLE
%       movePossible = E042.qVMO({'X'}, 10)
%       movePossible = E042.qVMO({'X' 'Y' 'Z'}, [10 2.2 0.8])
%       movePossible = E042.qVMO(E042.qSAI, [-10 -10 -10 10 10 10 10 10])
%
%   E042 is an alias for an arbitrary controller name.

% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(c.ID<0),                             error('The controller is not connected'), end;
if (length(axes) ~= length(positions)), error('values vector size does not fit identifier string size'), end;
szAxes = CellArrayToStringSpaceSeperated(axes);
functionName = [ c.libalias, '_', mfilename];

if(any(strcmp(functionName,c.dllfunctions)))
	pdValues = libpointer('doublePtr', positions);
    movePossible = zeros(1,1);
	piMovePossible = libpointer('int32Ptr',movePossible);
	try
		[bRet, ~, ~, movePossible] = calllib(c.libalias, functionName, c.ID, szAxes, pdValues, piMovePossible);
		if(bRet==0)
			iError = GetError(c);
			szDesc = TranslateError(c,iError);
			error(szDesc);
		end
    catch ME
        error (ME.message);
	end
else
	error('%s not found',functionName);
end

movePossible = movePossible==1;
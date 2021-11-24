function resPos = VFN_PIStage_move(stage, pos, timeout)
% VFN_PIStage_move Function to move a PI stage IN CLOSED LOOP
%   Values handled in [mm] 
%
% This will block execution until the reference completes
%   a user-provided timeout is implemented on the block
%   if user doesn't provide a timeout, default of 120s is used
%   The timeout throws an easy-to-catch error for error handling
%
% NOTE: assumes the axis name is '1'
% NOTE: assumes servo is already on:
%       - use VFN_PIStage_setServo() to turn on if needed
% NOTE: assumes the axis has already been referenced:
%       - use VFN_PIStage_home() to reference if needed
%   
%   EXAMPLE:______
%   resPos = VFN_PIStage_move(stage, pos)
%       'stage'    = instance of PI device (USB connection) object
%                    (ie. a single field in the PIdevs struct)
%       'pos'      = desired position in range (-13,13) in [mm]
%       'timeout'  = (OPTIONAL) time to wait [sec] for move to complete
%                    defualt timeout = 120
%       'resPos'   = resulting position in [mm]


                        
% Assume axis name is '1'
axis = '1';

% Provide default value for timeout if needed
if nargin < 3
    timeout = 120;
end
                        
%% Perform move in closed loop
stage.MOV(axis, pos);

%% Wait for move to finish
tic     % implement timeout to avoid infinite loop
while(stage.IsMoving == 1)
    pause(0.01);
    if toc > timeout
        error('VFN:PIStage:move_timeout','Stage failed to finish move in %d seconds',timeout)
    end
end 

%% Get resulting position
resPos = stage.qPOS(axis);

end
function resPos = VFN_PIStage_move(stage, pos)
% VFN_PIStage_move Function to move a PI stage IN CLOSED LOOP
%   Values handled in [mm] 
%
% This will block execution until the reference completes
%   a 60s timeout is implemented on the block
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
%       'resPos'   = resulting position in [mm]


                        
% Assume axis name is '1'
axis = '1';
                        
%% Perform move in closed loop
stage.MOV(axis, pos);

%% Wait for move to finish
tic     % implement timeout to avoid infinite loop
while(stage.IsMoving == 1)
    pause(0.01);
    if toc > 60
        error('Stage failed to finish move in 60 seconds')
    end
end 

%% Get resulting position
resPos = stage.qPOS(axis);

end
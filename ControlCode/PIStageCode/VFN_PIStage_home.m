function VFN_PIStage_home(stage)
% VFN_PIStage_home Function to home a PI stage
% 
% This will set the servo switch on the stage to ON
% This will block execution until the reference completes
%   a 60s timeout is implemented on the block
%
% NOTE: assumes the axis name is '1'
%   
%   EXAMPLE:______
%   VFN_PIStage_home(stage)
%       'stage' = instance of PI device (USB connection) object
%                 (ie. a single field in the PIdevs struct)
%


% Assume axis name is '1'
axis = '1';

%% Confirm argument is valid
if ~isa(stage, 'PI_GCS_Controller')
    error("Argument must be a 'PI_GCS_Controller' object but was a '%s'", class(stage))
end

%% Switch servo on (needed for home move)
stage.SVO(axis, 1);

%% Reference the axis
stage.FRF(axis);  % perform reference

disp('Referencing (homing) stage')

%% Wait for referencing to finish
tic
while(stage.qFRF(axis) == 0)                        
    pause(0.5);
    if toc > 60
        error('Stage failed to home in 60 seconds')
    end
end

end

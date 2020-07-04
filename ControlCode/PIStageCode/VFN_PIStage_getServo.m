function state = VFN_PIStage_getServo(stage)
% VFN_PIStage_getServo Function to query the servo state on a PI stage
%   When the "servo" is on, the stage is in Closed Loop.
%   When the "servo" is off, the stage is in Open Loop.
% 
% NOTE: assumes the axis name is '1'
%   
%   EXAMPLE:______
%   state = VFN_PIStage_getServo(stage)
%       'stage' = instance of PI device (USB connection) object
%                 (ie. a single field in the PIdevs struct)
%       'state' = (number) 1='on'; 0='off'
%


% Assume axis name is '1'
axis = '1';

%% Confirm argument is valid
if ~isa(stage, 'PI_GCS_Controller')
    error("Argument must be a 'PI_GCS_Controller' object but was a '%s'", class(stage))
end

%% Set the servo state
state = stage.qSVO(axis);

end
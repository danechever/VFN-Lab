function VFN_PIStage_setServo(stage, state)
% VFN_PIStage_setServo Function to set the servo on a PI stage on/off
% 
% NOTE: assumes the axis name is '1'
%   
%   EXAMPLE:______
%   VFN_PIStage_setServo(stage, state)
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
stage.SVO(axis, state);

end
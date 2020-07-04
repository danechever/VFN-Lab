function res = VFN_PIStage_isHomed(stage)
% VFN_PIStage_isHomed Function to check if a PI stage is referenced/homed
% 
% NOTE: assumes the axis name is '1'
%   
%   EXAMPLE:______
%   res = VFN_PIStage_isHomed(stage)
%       'stage' = instance of PI device (USB connection) object
%                 (ie. a single field in the PIdevs struct)
%       'res'   = 1=reference valid; 0=not referenced
%


% Assume axis name is '1'
axis = '1';

%% Confirm argument is valid
if ~isa(stage, 'PI_GCS_Controller')
    error("Argument must be a 'PI_GCS_Controller' object but was a '%s'", class(stage))
end

%% Check if stage is referenced
res = stage.qFRF(axis);

end

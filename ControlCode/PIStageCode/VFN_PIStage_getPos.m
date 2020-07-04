function pos = VFN_PIStage_getPos(stage)
% VFN_PIStage_getPos Function for quering the position of a PI stage
%   Value returned in [mm]
%
% NOTE: assumes the axis name is '1'
%   
%   EXAMPLE:______
%   pos = VFN_PIStage_getPos(stage)
%       'stage' = instance of PI device (USB connection) object
%                 (ie. a single field in the PIdevs struct)
%       'pos'   = output position in [mm]


% NOTE: I removed lots of error checking to reduce runtime

% Assume axis name is '1'
axis = '1';

%% Get position
pos = stage.qPOS(axis);

end
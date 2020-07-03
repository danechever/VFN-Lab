function resPos = VFN_Zab_getPos(axis)
% VFN_Zab_resPos Function for quering the position of a Zaber axis
%   
%   EXAMPLE:______
%   resPos = VFN_Zab_getPos(axis)
%       axis:       Axis to home. Must be an instance of zaber.motion.ascii.Axis
%       resPos:     resulting position (in mm)

if ~isa(axis, 'zaber.motion.ascii.Axis')
    error('The first argument must be a zaber.motion.ascii.Axis object')
end

%% get the position of each axis
% Query position; convert to mm
try
    resPos = axis.getPosition(zaber.motion.Units.LENGTH_MILLIMETRES);
catch exception
    % Close port if a MATLAB error occurs, otherwise it remains locked
    axis.getDevice().getConnection().close();
    rethrow(exception);
end


end
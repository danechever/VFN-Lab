function resPos = VFN_Zab_moveRel(axis, step)
% VFN_Zab_moveRel Function for performing a relative move with a Zaber axis
%   Blocks execution until move completes
%   
%   EXAMPLE:______
%   resPos = VFN_Zab_moveRel(axis, pos)
%       axis:       Axis to home. Must be an instance of zaber.motion.ascii.Axis
%       step:       Size of relaitve move (in mm)
%       resPos:     resulting position (in mm)

if ~isa(axis, 'zaber.motion.ascii.Axis')
    error('The first argument must be a zaber.motion.ascii.Axis object')
end

%% Move axis
if ~isnan(step)
    % Perform move and capture device error code
    try
        axis.moveRelative(step, zaber.motion.Units.LENGTH_MILLIMETRES);
    catch exception
        % Close port if a MATLAB error occurs, otherwise it remains locked
        axis.getDevice().getConnection().close();
        rethrow(exception);
    end
end

%% Get Final Position

% Query position; convert to mm
try
    resPos = axis.getPosition(zaber.motion.Units.LENGTH_MILLIMETRES);
catch exception
    % Close port if a MATLAB error occurs, otherwise it remains locked
    axis.getDevice().getConnection().close();
    rethrow(exception);
end

end
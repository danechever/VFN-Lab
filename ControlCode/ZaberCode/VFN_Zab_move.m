function resPos = VFN_Zab_move(axis, pos)
% VFN_Zab_move Function for moving a Zaber axis
%   Blocks execution until move completes
%
% NOTE: If an out-of-range move is requested, the port will be closed and
% you will need to reopen it.
%
%   EXAMPLE:______
%   resPos = VFN_Zab_move(axis, pos)
%       axis:       Axis to move. Must be an instance of zaber.motion.ascii.Axis
%       pos:        position to move to (in mm)
%       resPos:     resulting position (in mm)

if ~isa(axis, 'zaber.motion.ascii.Axis')
    error('The first argument must be a zaber.motion.ascii.Axis object')
end

%% Move axis
if ~isnan(pos)
    % Perform move and capture device error code
    try
        axis.moveAbsolute(pos, zaber.motion.Units.LENGTH_MILLIMETRES);
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
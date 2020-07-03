function VFN_Zab_home(axis)
% VFN_Zab_resPos Function for homing a Zaber axis
%   
%   EXAMPLE:______
%   resPos = VFN_Zab_move(axis)
%       axis:       Axis to home. Must be an instance of zaber.motion.ascii.Axis

if ~isa(axis, 'zaber.motion.ascii.Axis')
    error('The first argument must be a zaber.motion.ascii.Axis object')
end

try
    axis.home();
catch exception
    % Close port if a MATLAB error occurs, otherwise it remains locked
    axis.getDevice().getConnection().close();
    rethrow(exception);
end

end
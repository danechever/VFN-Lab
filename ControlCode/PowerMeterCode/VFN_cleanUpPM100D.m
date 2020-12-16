% Close the connection to the PM100D
%
%   - Also deletes PM variable from workspace
%
%   See also VFN_setUpPM100D, VFN_PM100D_read


if isfield(PM, 'rm')
    % Close connection to the device
    PM.rm.close();
end
% NOTE: no close needed for the usbtmc connectionc

% Delete PM variable
clear('PM')
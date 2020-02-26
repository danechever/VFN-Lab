% Close the connection to the PM100D
%
%   - Also deletes PM variable from workspace
%
%   See also VFN_setUpPM100D, VFN_PM100D_read


% Close connection to the device
PM.rm.close();

% Delete PM variable
clear('PM')
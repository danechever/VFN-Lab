% Close the connection to two simultaneous PM100Ds
%
%   - Also deletes PM variables from workspace
%
%   See also VFN_setUpTWOPM100Ds, VFN_PM100D_read

%% Close connections
if isfield(normPM, 'rm')
    % Close connection to the device
    normPM.rm.close();
end
if isfield(nullPM, 'rm')
    % Close connection to the device
    nullPM.rm.close();
end
% NOTE: no close needed for the usbtmc connectionc

%% Cleanup workspace
% Delete PM variable
clear('normPM','nullPM')
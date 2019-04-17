function VFN_Zab_home(axis)
% VFN_Zab_resPos Function for homing a Zaber axis
%   
%   EXAMPLE:______
%   resPos = VFN_Zab_move(axis)
%       axis:       Axis to home. Must be an instance of Zaber.AsciiDevice.

if ~isa(axis, 'Zaber.AsciiDevice')
    error('The first argument must be a Zaber Ascii Device object')
end

try
    axis.home();
catch exception
    % Close port if a MATLAB error occurs, otherwise it remains locked
    fclose(axis.Protocol.Port);
    rethrow(exception);
end

axis.waitforidle();

end
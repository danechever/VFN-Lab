function resPos = VFN_Zab_getPos(axis)
% VFN_Zab_resPos Function for quering the position of a Zaber axis
%   
%   EXAMPLE:______
%   resPos = VFN_Zab_getPos(axis)
%       axis:       Axis to query. Must be an instance of Zaber.AsciiDevice
%       resPos:     resulting position (in mm)

if ~isa(axis, 'Zaber.AsciiDevice')
    error('The first argument must be a Zaber Ascii Device object')
end

M2MM = 1000;    %Conversion by which to multiply m to get mm

%% get the position of each axis
% Query position; convert to mm
try
    resPos = axis.Units.nativetoposition(axis.getposition)*M2MM;
catch exception
    % Close port if a MATLAB error occurs, otherwise it remains locked
    fclose(axis.Protocol.Port);
    rethrow(exception);
end


end
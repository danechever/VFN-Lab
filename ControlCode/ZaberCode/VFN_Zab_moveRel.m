function resPos = VFN_Zab_moveRel(axis, step)
% VFN_Zab_moveRel Function for performing a relative move with a Zaber axis
%   
%   EXAMPLE:______
%   resPos = VFN_Zab_moveRel(axis, pos)
%       axis:       Axis to move. Must be an instance of Zaber.AsciiDevice.
%       step:       Size of relaitve move (in mm)
%       resPos:     resulting position (in mm)

if ~isa(axis, 'Zaber.AsciiDevice')
    error('The first argument must be a Zaber Ascii Device object')
end

M2MM = 1000;    %Conversion by which to multiply m to get mm

%% Move axis
% Convert step size from mm to m (required for Zaber library)
step = step/M2MM;

if ~isnan(step)
    % Perform move and capture device error code
    try
        err = axis.moverelative(axis.Units.positiontonative(step));
    catch exception
        % Close port if a MATLAB error occurs, otherwise it remains locked
        fclose(axis.Protocol.Port);
        rethrow(exception);
    end
    % Throw error if a device error occurred
   	if ~isempty(err)
        error('An error occurred while moving the zaber\n  Error code = %d', err)
    end
end

%% Get Final Position
%Wait for axis to finish moving
axis.waitforidle();

% Query position for output; convert to mm
try
    resPos = axis.Units.nativetoposition(axis.getposition)*M2MM;
catch exception
    resPos = axis.Units.nativetopos
    % Close port if a MATLAB error occurs, otherwise it remains locked
    fclose(axis.Protocol.Port);
    rethrow(exception);
end

end
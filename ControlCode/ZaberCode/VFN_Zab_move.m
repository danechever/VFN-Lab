function resPos = VFN_Zab_move(axis, pos)
% VFN_Zab_move Function for moving a Zaber axis
%   Blocks execution until move completes
%
%   EXAMPLE:______
%   resPos = VFN_Zab_move(axis, pos)
%       axis:       Axis to move. Must be an instance of Zaber.AsciiDevice.
%       pos:        position to move to (in mm)
%       resPos:     resulting position (in mm)

if ~isa(axis, 'Zaber.AsciiDevice')
    error('The first argument must be a Zaber Ascii Device object')
end

M2MM = 1000;    %Conversion by which to multiply m to get mm

%% Move axis
% Convert position values from mm to m (required for Zaber library)
pos = pos/M2MM;

% Move Vertical axis if needed
if ~isnan(pos)
    % Perform move and capture device error code
    try
        err = axis.moveabsolute(axis.Units.positiontonative(pos));
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
    % Close port if a MATLAB error occurs, otherwise it remains locked
    fclose(axis.Protocol.Port);
    rethrow(exception);
end

end
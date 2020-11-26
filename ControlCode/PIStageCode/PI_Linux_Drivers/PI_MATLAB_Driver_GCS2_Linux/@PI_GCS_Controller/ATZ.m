function ATZ(c, szAxesIdentifierArray, lowVoltageArray)
% AuTo Zero: Automatic zero-point calibration.
%
%   SYNTAX
%       devices = ATZ(szAxes)
%       devices = ATZ(szAxes, lowVoltageArray)
%
%   DESCRIPTION
%   AuTo Zero: Automatic zero-point calibration. Sets the value of a
%   property (e.g. output voltage or axis position) which is to be applied
%   when the value of the axis property to be calibrated is zero. Starts an
%   appropriate calibration procedure.
%   The axis property to be calibrated with ATZ depends on the controller
%   and can be, for example, position or force.
%   This command can be interrupted by STP.
%   ATZ works in open-loop operation (servo off). If the servo is on, it
%   will be switched off automatically at the start of the ATZ procedure
%   and switched on again when it is finished.
%   The AutoZero procedure has the highest priority, i.e. it will overwrite
%   the control values given by all other sources. When the analog control
%   input is enabled, it will be disabled automatically at the start of the
%   AutoZero procedure and reenabled again when AutoZero is finished.
%   ATZ is not effective on non-linear axes (rotation axes) and walking
%   drives.
%   The success of the automatic zero-point calibration can be queried with
%   the qATZ command. Note that starting the AutoZero procedure for an
%   individual axis may influence the AutoZero results of other axes so
%   that their success states are reset. For that reason it is recommended
%   to start the AutoZero procedure for all axes at the same time (see
%   below).
%   The automatic zero-point calibration can take several seconds. During
%   this time, the controller is busy and only very limited able to execute
%   or answer commands. The ATZ procedure will move the axis, and the
%   motion may cover the whole travel range. Make sure that it is safe for
%   the stage to move.
%
%   INPUT
%       axesIdentifierArray   	String with axes
%
%       pdLowVoltageArray       Array with low voltages for the corresponding
%                               axes. If this parameter is omitted, the value
%                               stored in the controller (Autozero Low Voltage
%                               parameter, ID 0x07000A00) is used.
%
%   EXAMPLE
%       devices = ATZ()         #todo
%       devices = ATZ('E-042')  #todo
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Settings


%% Function arguments handle

numberOfAxes = GetNrAxesInString(c, szAxesIdentifierArray);

if nargin<3
    useDefaultArray = ones(numberOfAxes, 1);
    lowVoltageArray = zeros(numberOfAxes, 1);
else
    useDefaultArray = zeros(numberOfAxes, 1);
    
    if (numberOfAxes ~= length(lowVoltageArray))
        error('The two input arguments do not have the same length. "axesIdentifierArray" and "lowVoltageArray" must have the same length or "lowVoltageArray" must be omitted in (see help for ATZ).');
    end
end


%% Check input parameters


%%  Program start

functionName = [ c.libalias, '_', mfilename];

if (any(strcmp(functionName,c.dllfunctions)))
    pdLowVoltageArray = libpointer('doublePtr',lowVoltageArray);
    pbUseDefaultArray = libpointer('int32Ptr',useDefaultArray);
    try
        [bRet] = calllib(c.libalias,functionName, c.ID, ...
            szAxesIdentifierArray,pdLowVoltageArray,pbUseDefaultArray);
        if (bRet == 0)
            iError = GetError(c);
            szDesc = TranslateError(c,iError);
            error(szDesc);
        end
    catch
        rethrow(lasterror);
    end
else
    error('%s not found',functionName);
end
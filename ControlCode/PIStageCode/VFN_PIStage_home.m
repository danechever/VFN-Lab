%% Comment:
% This will home the "axis" on an E-873 controller
% It assumes the following variables have already been created by _init.m
    % 1) 'axis'         = char containing the axis name for use by other scripts
    % 2) 'PIdevice'     = instance of device (USB connection) object
% It will set the the servo switch on the stage to "ON"
% This command will block until the reference completes
    % NOTE: a 60s timeout is implemented on the block

% switch servo on for axis
switchOn    = 1;
% switchOff   = 0;
PIdevice.SVO ( axis, switchOn );

% reference axis
PIdevice.FRF ( axis );  % find reference

fprintf('Referencing (homing) axis %s\n',axis)
% wait for referencing to finish
tic     % implement timeout to avoid infinite loop
ncnt    = 1;    %counter for time display
fprintf('-- time elaped [s]: 0; ')
while(0 ~= PIdevice.qFRF ( axis ) == 0 )                        
    pause(0.1);           
    %[char ( 8 ) * ones( 1, 7 ), '.']
    if ~mod(ncnt, 50)
        % approx. 5s have passes, print so
        fprintf('%.0f; ', ncnt/10)
    end
    ncnt = ncnt + 1;
    if toc > 60
        error('\nStage on axis %s failed to home in 60 seconds',axis)
    end
end 
% print new line to clear time display line
fprintf('\nAxis %s referenced (homed)\n', axis)

clear switchOn ncnt
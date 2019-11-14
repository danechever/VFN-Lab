%% Comment:
% This will move the "axis" on an E-873 controller to position "PIpos" 
% It assumes the following variables have already been created by _init.m
    % 1) 'axis'         = char containing the axis name for use by other scripts
    % 2) 'PIdevice'     = instance of device (USB connection) object
% It also assumes the following variable is defined elsewhere
    % 3) 'PIpos'        = position to which the device should be moved
                        %   This should be a float in range (-13, 13) [mm]
% This command will block until the move completes
    % NOTE: a 60s timeout is implemented on the block

%% Perform move
%-- Check that axis has been referenced
if ~PIdevice.qFRF(axis)
    error('Axis %s has not been referenced (home). Call _home.m', axis)
end

%-- Check that servo is on for closed-loop move
svoFlag = false;
if ~PIdevice.qSVO(axis)
    svoFlag = true;
    PIdevice.SVO ( axis, 1);
    warning('Servo was off for axis %s. Servo has been turned on for move',axis)
end

%-- Move in closed loop
PIdevice.MOV ( axis, PIpos);

% OLD BLOCKING FEATURE -- too verbose and inneficient
% disp ( 'Stage is moving')
% % wait for motion to stop
% while(0 ~= PIdevice.IsMoving ( axis ) )
%     pause ( 0.1 );
%     [char ( 8 ) * ones( 1, 7 ), '.']
% end

% NEW BLOCKING FEATURE -- longer but more efficient
fprintf('Moving axis %s\n',axis)
% wait for move to finish
tic     % implement timeout to avoid infinite loop
ncnt    = 1;    %counter for time display
while(0 ~= PIdevice.IsMoving ( axis ) )                        
    pause(0.01);           
    %[char ( 8 ) * ones( 1, 7 ), '.']
    if ~mod(ncnt, 500)
        if ncnt == 500
            %Start printing on first call if needed
            fprintf('-- time elaped [s]: ')
        end
        % approx. 5s have passes, print so
        fprintf('%.0f; ', ncnt/100)
    end
    ncnt = ncnt + 1;
    if toc > 60
        error('\nStage on axis %s failed to home in 60 seconds',axis)
    end
end 
% print new line to clear time display line
fprintf('\n')

%-- Print final position
finPos = PIdevice.qPOS(axis);
fprintf('Axis %s moved to %f\n', axis, finPos)

if svoFlag
    PIdevice.SVO ( axis, 0);
    warning('Servo on axis %s was turned back off since it was off before the move', axis)
end

clear svoFlag ncnt finPos
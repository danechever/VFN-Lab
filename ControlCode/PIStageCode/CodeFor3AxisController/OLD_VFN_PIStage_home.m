function VFN_PIStage_home(axis, PIdevice)
%% Comment:
% This will home the "axis" on an E-873 controller
% It assumes the following variables have already been created by _init.m
    % 1) 'axis'         = mutlidimensional cell containing the axis names for use by other scripts
    % 2) 'PIdevice'     = instance of device (USB connection) object
% It will set the the servo switch on the stage to "ON"
% This command will block until the reference completes
    % NOTE: a 60s timeout is implemented on the block

% switch servo on for axis
switchOn    = ones(1,length(axis));
% switchOff   = 0;
PIdevice.SVO ( axis, switchOn );


% reference axis
PIdevice.FRF ( axis );  % find reference
ref_axis = containers.Map;
ref_axis('1') = 'x'; ref_axis('2') = 'y'; ref_axis('3') = 'z';
% axis 1 is x; axis 2 is y; axis 3 is z

for axisID=1:length(axis)
    fprintf('Referencing (homing) axis %s\n', char(ref_axis(char(axis(axisID)))))
end

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
        fprintf('Most recent call involved stages:\n')
        for axisID = 1:length(axis)
            fprintf('%s\n', ref_axis(char(axis(axisID))))
        end
        error('\nStage in most recent call failed to home in 60 seconds')
    end
end 
% print new line to clear time display line
for axisID=1:length(axis)
    fprintf('\nAxis %s referenced (homed)\n', char(ref_axis(char(axis(axisID)))))
end

clear switchOn ncnt ref_axis
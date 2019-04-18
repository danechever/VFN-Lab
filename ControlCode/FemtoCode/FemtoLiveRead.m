isAutoScale = false; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO_scale = 5;     % Starting gain. the n in: 10^n for the gain setting
Nread   = 100;      % power samples at a given locaiton
Nrate   = 1000;     % rate of scan [samples/second]
Delay = 0.25;

%% Femto setup
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

%if isAutoScale
%    % Take a sample reading at current power
%    readVal = mean(startForeground(s));
%    % Modify gain accordingly
%    [FMTO_scale, s] = VFN_FMTO_setAutoGain(s, readVal, FMTO_scale);
%else
    s = VFN_FMTO_setGain(s, FMTO_scale);
%end

fprintf('Delay in use: %0.1f\n', Delay)
fprintf('Current Gain setting: %i\n', FMTO_scale)

%% Read Femto
for i = 1:10000
    read    = startForeground(s);
    if isAutoScale
        % Save old FMTO_scale for comparison after autoGain
        old_scale = FMTO_scale;
        % Save old read value for comparison after autoGain
        old_read = mean(read);
        % Modify gain accordingly
        [FMTO_scale, s] = VFN_FMTO_setAutoGain(s, old_read, FMTO_scale);
        % Check if re-read is needed
        if old_scale ~= FMTO_scale
            % FMTO_scale changed so the gain changed
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
            read    = startForeground(s);
            %ratio   = ratio * (old_read/mean(read));
            if FMTO_scale > 9
                warning('Gain >9: %i',FMTO_scale)
            end
        end
    end
    if find(read>10)
        warning('Power is too high')
    end
    fprintf('Itr %05d Gain: %d Power: %f\n', i, FMTO_scale, mean(read))
    pause(Delay)
end
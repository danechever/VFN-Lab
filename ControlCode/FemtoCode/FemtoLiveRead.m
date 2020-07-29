FMTO.isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO.FMTO_scale = 5;     % Starting gain. the n in: 10^n for the gain setting
FMTO.Nread   = 100;      % power samples at a given locaiton
FMTO.Nrate   = 1000;     % rate of scan [samples/second]
Delay = 0.25;
totread = 50;    % Total number of samples to take
isplot  = true;     % Flag to plot full vector at end of run

%% Femto setup
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
FMTO.gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

VFN_FMTO_LUCI_setGain(FMTO);

fprintf('Delay in use: %0.1f\n', Delay)
fprintf('Current Gain setting: %i\n', FMTO.FMTO_scale)

if isplot
    meas = nan(totread,1);
end

%% Read Femto
for i = 1:totread
    read    = startForeground(FMTO.s);
    if FMTO.isAutoScale
        % Save old FMTO_scale for comparison after autoGain
        old_scale = FMTO.FMTO_scale;
        % Save old read value for comparison after autoGain
        old_read = mean(read);
        % Modify gain accordingly
        FMTO = VFN_FMTO_LUCI_setAutoGain(FMTO, old_read);
        % Check if re-read is needed
        if old_scale ~= FMTO.FMTO_scale
            % FMTO_scale changed so the gain changed
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO.FMTO_scale)
            read    = startForeground(FMTO.s);
            %ratio   = ratio * (old_read/mean(read));
            if FMTO.FMTO_scale > 9
                warning('Gain >9: %i',FMTO.FMTO_scale)
            end
        end
    end
    if find(read>10)
        warning('Power is too high')
    end
    fprintf('Itr %05d Gain: %d Power: %f\n', i, FMTO.FMTO_scale, mean(read))
    if isplot
        meas(i) = mean(read)*FMTO.gnFact^-(FMTO.FMTO_scale-6);
    end
    pause(Delay)
end

VFN_cleanUpFMTO;

%% Plot results if needed
if isplot
    figure; plot(meas);
end
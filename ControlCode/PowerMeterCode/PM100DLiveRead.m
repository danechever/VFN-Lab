timerd  = 16*60;       % [min] Time for which reads should be done (approx)
Nrate   = 50;       % number of measurements to average for each sample
Delay   = 60;       % [s] Amount of time between samples
isplot  = true;     % Flag to plot full vector at end of run

%% PM100D Setup

% Setup PM100D 
VFN_setUpPM100D;  


%% Read PM100D
%-- Calculate number of reads to perform
min2sec = 60;   % [s/min] number of seconds per minute
nRead = timerd*min2sec*(1/Delay);   % Total number of samples to take


fprintf('Total time to run:    %0.2f [min]\n', timerd)
fprintf('Time between samples: %0.2f [s]\n', Delay)
fprintf('Number of samples:    %i\n', nRead)

if isplot
    meas = nan(nRead,Nrate);
end

for i = 1:nRead
    % Take measurements
    read    = VFN_PM100D_read(PM,Nrate);
    %read    = mean(read);
    
    % Print values
    fprintf('Itr %05d | Power: %f [mW]\n', i, mean(read)*1e3);
    if isplot
        % Save data for plot
        meas(i,:) = read;
    end
    pause(Delay)
end

VFN_cleanUpPM100D;

%% Plot results if needed
if isplot
    figure; plot(meas);
end

%-- OPTION 1: Run slowlish, sampling (roughly) at a specified time interval
% isOpt   = 1;        % Flag denoting option choice (do not change, just uncomment)
% timerd  = 0.5;       % [min] Time for which reads should be done (approx)
% Nrate   = 50;       % number of measurements to average for each sample
% Delay   = 3;       % [s] Amount of time between samples
% isplot  = true;     % Flag to plot full vector at end of run

%-- OPTION 2: Run as fast as possible, sampling whenever possible
isOpt   = 2;        % Flag denoting option choice (do not change, just uncomment)
timerd  = 2;        % [min] Time for which reads should be done (approx)
Nrate   = 50;       % number of measrements to average in each sample
Delay   = 0.3;     % [s] Measured minimum time between samples
isplot  = true;     % Flag to plot full vector at end of run

%% PM100D Setup

% Setup PM100D 
VFN_setUpPM100D;  


%% Read PM100D
%-- Calculate number of reads to perform
min2sec = 60;   % [s/min] number of seconds per minute
nRead = floor(timerd*min2sec*(1/Delay));   % Total number of samples to take


fprintf('Sampling choice:      %i\n', isOpt)
fprintf('Total time to run:    %0.2f [min]\n', timerd)
fprintf('Time between samples: %0.2f [s]\n', Delay)
fprintf('Number of samples:    %i\n', nRead)

if isplot
    meas = nan(nRead,Nrate);
end

% Preallocate time vector
timSamps = nan(nRead,1);

tic     % Start "timer" for time vector
for i = 1:nRead
    % Take measurements
    read    = VFN_PM100D_read(PM,Nrate);
    %read    = mean(read);
    
    % Print values
    fprintf('Itr %05d | Power: %0.8f [mW]\n', i, mean(read)*1e3);
    timSamps(i) = toc;
    if isplot
        % Save data for plot
        meas(i,:) = read;
    end
    if isOpt == 1
        pause(Delay)
    end
end
toc

VFN_cleanUpPM100D;

% Save results
%save('/media/Data_Drive/VFN/TestbedData/210428_COV8/superKCompact_650x10nm_pmRead_12hr_fastSamp_80percentPower.mat','timSamps','meas')

%% Plot results if needed
if isplot
    figure; plot(timSamps, mean(meas,2));
    drawnow; pause(0.05);
end
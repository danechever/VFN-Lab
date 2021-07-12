%{
Script for obtaining Femto Gain impulse response, rise time, and bias.

This script measures how long it takes for the femto power to settle after
changing the gain setting. It samples the femto immediately after a gain
change then displays the resulting voltage-vs-time curve. From that curve,
it computes the rise time required for the voltage to stabilize at a
steady-state value.

The bias at the given gain is reported as the steady-state value.

The script procedure is:
- Set desired femto gain setting to sample
- Measure power on femto immediately after and as quickly as possible
- Compute the steady-state voltage as the average from the later samples.
- Report the "bias" at the given gain as this steady-state voltage.
- Compute the rise time as the time required to go from the initial voltage
  to 95% of the steady-state voltage.

Optical setups should be:
- femto input port covered with the black rubber cap (so that it is
  measuring no light power).

There is no saved output. The script displays the results directly in the
command window and on plots.
%}

close all; clear all;
FMTO.FMTO_scale = 5;     % Starting gain. the n in: 10^n for the gain setting
FMTO.Nread   = 100;      % power samples in a single burst
Delay = 0.0;        % Time to wait between sample bursts
totread = 150;    % Total number of sample bursts to take
isplot  = true;     % Flag to plot full vector at end of run

RT_frac = 0.95;  % Fraction of steady-state value to consider for rise time
%% Femto setup
% Setup Femto 
VFN_setUpFMTO;  

fprintf('Delay in use: %0.1f\n', Delay)
fprintf('Current Gain setting: %i\n', FMTO.FMTO_scale)
fprintf('Sample Bursts to read: %d\n',totread)

% Preallocate data matrices
meas = nan(totread,FMTO.Nread);
times = nan(totread,1);

% Set isBlock = false so that we can see the voltage readings immediately
% after the gain change. (This is an impulse response script afterall)
VFN_FMTO_setGain(FMTO,false);
% Save Gain used
gainset = FMTO.FMTO_scale;
%pause(0)

% Start timmer
tic

%% Read Femto
for iii = 1:totread
    meas(iii,:) = VFN_FMTO_readN(FMTO);
    times(iii) = toc;
    %pause(Delay)
end
VFN_cleanUpFMTO;

fprintf('Done measuring\n')
disp('')
fprintf('Total Time:          %f\n',times(end))
fprintf('Avg Time per sample: %f\n', mean(diff(times)))
%% Get statistics
%-- Compute the mean and STD
meas_mean = mean(meas,2);
meas_std  = std(meas,[],2);

%-- Compute the steady-state value
  % Assume second half of measurements are steady-state
ss_mean = mean(meas_mean(end/2:end));
ss_std  = mean(meas_std(end/2:end));

%-- Compute the rise time
% First, offset all voltages so that the min value = 0
tmp_mean = meas_mean-min(meas_mean);
% Now find the index of value that crosses the rise time fraction provided
    % by the user
RT_ind = find(tmp_mean > (RT_frac*(ss_mean-min(meas_mean))),1);
% Get the rise resulting rise time
RT = times(RT_ind);

fprintf('Steady-State Mean:   %f\n',ss_mean)
fprintf('Steady-State STD:    %f\n',ss_std)
fprintf('Rise Time:           %f\n\n',RT)

%% Plot results if desired
if isplot
    figure(); 
    plot(times, meas_mean);
    xlabel('Time [s]')
    ylabel('Mean of Power [V]')
    title(sprintf('Power - Gain %d - SS Mean %f - RiseTime %f',gainset,ss_mean,RT))
    % add the ss value
    hold on
    yline(ss_mean)
    % add the risetime marker
    xline(RT)
    
    figure();
    plot(times, meas_std);
    xlabel('Time [s]')
    ylabel('STD of Power [V]')
    title(sprintf('Standard deviations - Gain %d',gainset))
    % add the ss value
    hold on
    yline(ss_std)
    % add the risetime marker
    xline(RT)
    
    figure();
    plot(times);
    xlabel('Sample')
    ylabel('Time [s]')
    title(sprintf('Time Sampling - Gain %d',gainset))
end
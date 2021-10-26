%{
Script for obtaining Femto-to-redPM Ratio - (Manually)

*** NOTE: THIS SCRIPT IS MEANT TO BE USED MANUALLY. It'll prompt you to set
    a new test condition and then will measure the power on both the femto
    and redPM in series.

*** This is meant to reduce the number of times the varia parameters are
    changed since it seems the varia may not be super repeatable. The old
    code (Femto_redPM_Bandwidth_Calibrator) would have to scan through all
    bandwidths twice: once to measure on the redPM and once to measure on
    the femto. This new code instead prompts the user to set a bandwidth on
    the varia manually before measuring the power on the redPM and then
    prompting the user to change the fiber to the femto.

*** This code could hypothetically be used to calibrate the redPM and femto
    on other light sources as well.

This script uses the femto and redPM only. The goal is to measure the femto
and redPM power at various bandwidths simultaneously. These measurements 
can be used to obtain the ratio between the femto voltage and redPM for a 
given bandwidth.

The script procedure is:
1. Connect to femto and redPM
2. Prompt user to set the light-source test-condition. Code tells the user
    what condition is expected (lower and upper bandedges)
3. Prompt user to connect the fiber into the Femto or redPM first
4. Measure the power on the Femto/redPM at the desired gain settings.
5. Prompt user to switch the fiber to the redPM/Femto
6. Measure the power on the redPM/Femto
7. Prompt user to set new light-source test-condition. Now repeat #3-7
    until all test conditions are sampled.

*** NOTE: We alternate between reading from the femto or redPM first at each
bandwidth so that it minimizes the number of times the fiber needs to be
switched.

Optical setups should be:
- Light source output plugged into either the redPM or the femto with both
    available within reach of the optical fiber.

The output is 2 matrices of power/voltage readings, one for the redPM and
one for the Femto, where each row represents a different BW setting and
each column represents a gain setting. NOTE: since the redPM has no gain
settings, redPM data matrices are just a column vector.

NOTE: the script avoids sampling gains that are/will be saturated.
%}


close all; clear all;
addpath(genpath(['..']));

%% Define test conditions

% Choose which measuremente device to sample first
measdev0 = 'femto';   % must be either "Femto" or "redPM"

% Define number of sample points in a given burst of data
Nread   = 100;

% Define the gain values to try. Each value is the n in 10^n gain settings
  % Values in this vector should be in assending order
  % THIS ONLY APPLIES TO THE FEMTO
gainVals = [5, 6, 7, 8, 9, 10, 11];     

% Enter the light source settings to use for the test
  % NOTE: the code won't set the varia to these conditions, this is just so
  % that it can prompt the user to use the right parameters.
LAScentL = 650;     % [nm] central wavelength (single value only)
% Enter the bandwidths
  % NOTE: MAKE SURE YOU GET THE RIGHT ORDER
LASBWs   = flip([6, 10:10:100]);  % [nm] Bandwidths to sample

% Define the time to wait at various points in the code (settling time for
% the femto gain, settling time on laser changes, etc.)
Delay = 2.5;

% Define the value to use when a datapoint should be ignored
wierdVal1 = -8888;   % Value for when the current gain setting is saturated
wierdVal2 = -9999;   % Value for when a previous gain setting was saturated

% Provide filename for saving data
svnm = sprintf('/media/Data_Drive/VFN/TestbedData/210428_COV8/Femto2redPMCalibration_onSuperKEXT_BWs_RedoForNewData_Sample5',measdev0);
% Flag marking whether to save data or not
isSave = true;

%% Connect to the redPM and Femto

%-- Make sure that the user provided a valid measdev0 string
if ~strcmpi(measdev0, 'redpm') && ~strcmpi(measdev0, 'femto')
    error('The provided measdev0 was not recognized')
end

%-- Femto setup
disp('Setting up femto')
% Set the number of reads to do in a single burst
FMTO.Nread   = Nread;      % power samples in a single burst

% Connect to the femto
VFN_setUpFMTO;

%-- redPM setup
disp('Setting up redPM')

% Connect to the redPM
VFN_setUpPM100D

% Set operating wavelength for calibration
VFN_PM100D_setWavelength(PM, LAScentL);

% Wait a moment just in case the wavelength setting needs to settle
pause(1)

%% Take Data
%-- Preallocate measurements data matrix
  % Axes: (BW setting, samples in a single burst, gainValues)
meas_femto = nan(length(LASBWs), Nread, length(gainVals));
meas_redpm = nan(length(LASBWs), Nread);

%-- Set current measurement device to the user-provided measdev0 for start
measdev = measdev0;

%-- Iterate through Bandwidth settings and take measurements
for BWind = 1:length(LASBWs)
    fprintf('---- Iteration %d of %d ----\n',BWind,length(LASBWs))
    %-- Prompt user to set wavelengths
    % Compute the settings the user should use
    varlow = LAScentL-LASBWs(BWind)/2;
    varhi  = LAScentL+LASBWs(BWind)/2;
    % Prompt user to set those settings    
    fprintf('*** Set light source to [%0.2f , %0.2f]\n   Press any key when ready\n',varlow,varhi);
    % Audible cue to notify user we're waiting for them
    beep;
    % Wait for user to press any key notifying that code should procede
    pause;    
    
    % Wait for superK to stabilize, just to be safe
    pause(Delay)
    
    %-- Use a for-loop to measure both devices in the right order. 
      % Not the cleanest way to do this but easy to implement
    for ind = 1:2
        if strcmpi(measdev, 'femto')
            %-- Prompt user to switch to femto
            fprintf('*** Switch fiber to Femto ***\n')
            % Audible cue to notify user we're waiting for them
            beep;
            % Wait for user to press any key notifying that code should procede
            pause;    
            
            %-- Deal with Femto measurement

            % Define a flag that lets us skip higher gain settings when a lower
              % setting has already saturated at 10V
            gainSatFlag = false;

            fprintf(' Try gain = ')
            for gind = 1:length(gainVals)
                fprintf('%d ', gainVals(gind))
                %-- Check that the previous gain setting wasn't saturated
                if gainSatFlag
                    % If gainFlag is true, a prior gain setting was already
                    % saturated so let's not bother measuring this gain setting

                    % Let's provide a specific value in the data matrix though 
                    % so that we know this setting was considered
                    meas_femto(BWind, :, gind) = wierdVal2;

                    fprintf('sat, ')

                    % Now skip the rest of this for-loop iteration 
                    continue
                end

                %-- Set the desired gain
                FMTO.FMTO_scale = gainVals(gind);
                % Set isBlock=false so that we can use the user-provided delay here
                VFN_FMTO_setGain(FMTO,false);
                % Wait the desired time for the value to settle, just to be safe
                pause(Delay)

                %-- Take Measurements
                % Read power on the Femto
                read = VFN_FMTO_readN(FMTO);

                % Check that the power is less than 10
                if find(read>10)
                    % The power was above 10. Handle accordingly

                    % First, set the gain setting down so that we don't leave the
                    % femto saturated for too long
                    FMTO.FMTO_scale = 5;
                    VFN_FMTO_setGain(FMTO);

                    % Provide a specific value in the data matrices so we know this
                    % setting was considered but found to be too large
                    meas_femto(BWind, :, gind) = wierdVal1;

                    fprintf('sat, ')

                    % Set the saturation flag true so that we skip any higher gains
                    gainSatFlag = true;
                else
                    % The power was below 10 V and we can keep running as normal
                    % Save the result
                    meas_femto(BWind, :, gind) = read;
                end
            end
            %-- Once we're done with the femto, change the measdev variable
            %  to redpm ONLY if this was the first sample at this bandwidth.
            %  If it was the second, we don't want to change it since that
            %  adds unnecessary fiber switches.
            if ind == 1
                measdev = 'redpm';
            end
            disp(' ')  % Clear the line in the command prompt
        elseif strcmpi(measdev,'redpm')
            %-- Prompt user to switch to femto
            fprintf('*** Switch fiber to redPM ***\n')
            % Audible cue to notify user we're waiting for them
            beep;
            % Wait for user to press any key notifying that code should procede
            pause;    
            
            %-- Deal with redPM measurement

            % Wait a moment before measurement just to be extra super safe
              % yes, we already waited. Oh well, wait more
            pause(Delay)

            % Measure the redPM 
            read = VFN_PM100D_read(PM, Nread);

            % Save the result
            meas_redpm(BWind, :)   = read;
            %-- Once we're done with the redPM, change the measdev variable
            %  to femto ONLY if this was the first sample at this bandwidth.
            %  If it was the second, we don't want to change it since that
            %  adds unnecessary fiber switches.
            if ind == 1
                measdev = 'femto';
            end
        end
    end
end

%% cleanup all connections
%-- Disconnect from measurement devices
% Set femto gain to lowest setting
FMTO.FMTO_scale = 5;
VFN_FMTO_setGain(FMTO);
% Close femto connection
VFN_cleanUpFMTO;
% Close redPM connection
VFN_cleanUpPM100D;

%% Do some analysis
%-- Check if any matrix elements remained unassigned
if any(isnan(meas_redpm(:)))
    error('Some element in redPM data matrix is nan')
end
if any(isnan(meas_femto(:)))
    error('Some element in femto data matrix is nan')
end

%-- Make a copy of the data for now
meas_redpm_raw = meas_redpm;
meas_femto_raw = meas_femto;

%-- Change any of the "wierdVal#' entries to nan
meas_redpm(meas_redpm == wierdVal1) = nan;
meas_redpm(meas_redpm == wierdVal2) = nan;
meas_femto(meas_femto == wierdVal1) = nan;
meas_femto(meas_femto == wierdVal2) = nan;

%-- Compute the STD along third dimension
meas_redpm_std = std(meas_redpm, [], 2);
meas_femto_std = std(meas_redpm, [], 2);

%-- Take averages along second dimension 
meas_redpm = mean(meas_redpm, 2);
meas_femto = mean(meas_femto, 2);

%-- Display the final data
meas_redpm
squeeze(meas_femto)
%% Save the data
if isSave
    fprintf('\n\nSAVING DATA TO: %s\n', svnm)
    save(svnm, 'gainVals', 'LASBWs', 'LAScentL', 'Delay', ...
        'wierdVal1', 'wierdVal2', 'meas_redpm', 'meas_redpm_raw', ...
        'meas_femto', 'meas_femto_raw');
end

%{
Script for obtaining Femto Gain Factor

NOTE: this script is nearly identical to FemtoGainCalibrator.m except this
script prompts the user to set the power manually instead of connecting to
the KLS635. Thus, user can calibrate the gain on any light source.

This script uses the femto and redPM, with a manually set power source, to
measure the femto response at various power settings and obtain the voltage
ratio between different gain settings.

The script procedure is: 
- Prompt user to manually set laser power 
- Set femto gain 
- Measure power on femto and redPM 
- Increase gain (iterating until all gain settings are sampled) 
- Prompt user to manually increase laser power
- Iterate until user stops script

Optical setups should be:
- Power source plugged into input port of fiber splitter
- femto plugged into "tap" output port of fiber splitter
- redPM plugged into "signal" output port of fiber splitter

The output is a pair of matrices where columns are different gain settings
and rows are different power values. One matrix contains the the femto
measurements and the other contains the redPM measurements. The output
results can be analyzed with the FemtoGainCalibrator_Analysis script.

NOTE: the script avoids sampling gainst that are/will be saturated.
%}


close all; clear all;

%% Define test conditions

% Define the gain values to try. Each value is the n in 10^n gain settings
  % Values in this vector should be in assending order
gainVals = [5, 6, 7, 8, 9, 10, 11];     

% Define number of laser power values to sample
NlaserSamps = 20;

% Define the time to wait after changing the power setting
Delay = 2.5;

% Define the value to use when a datapoint should be ignored
wierdVal1 = -8888;   % Value for when the current gain setting is saturated
wierdVal2 = -9999;   % Value for when a previous gain setting was saturated

% Provide filename for saving data
svnm = '/media/Data_Drive/VFN/TestbedData/210428_COV8/GainCalibration_onSuperK_Sample3';
% Flag marking whether to save data or not
isSave = true;

%% Femto setup
disp('Setting up femto')
% Set the number of reads to do in a single burst
FMTO.Nread   = 100;      % power samples in a single burst

% Connect to the femto
VFN_setUpFMTO;  

%% redPM setup
disp('Setting up redPM')
pmNread = 100;      % Number of samples to take
pmCalWvl = 650; % Wavelength for redPM for normalization

% Connect to the redPM
VFN_setUpPM100D

% Set operating wavelength for calibration
VFN_PM100D_setWavelength(PM, pmCalWvl);

% Wait a moment just in case the wavelength setting needs to settle
pause(1)

%% Take Data
%-- Preallocate data matrices
% Define matrix to hold femto measurements. 
   % Axes: (powerValue, gainValue, samples in a single burst)
FMTOmeas = nan(NlaserSamps, length(gainVals), FMTO.Nread);
% Define matrix to hold redPM measurements. 
   % Axes: (powerValue, femto gainValue, samples in a single burst)
PMmeas = nan(NlaserSamps, length(gainVals), pmNread);
% Define vector to hold laser settings
pwrVals = nan(NlaserSamps, 1);

%-- Iterate through power settings and take measurements
for pind = 1:NlaserSamps
    % Prompt the user to set the new power and then enter that value
      % Save the value into the laser power vector so we can reference it
      % later. Also add audible cue to help user stay aware of prompt.
    fprintf('---- Iteration %d of %d ----\n',pind,NlaserSamps)
    beep;
    pwrVals(pind) = input('Set Laser power then enter the value: ');
    
    % Wait for power to stabilize, just to be safe
    pause(2)
    
    fprintf('\n Laser power set to: %f %%', pwrVals(pind))
    
    % Define a flag that lets us skip higher gain settings when a lower
      % setting has already saturated at 10V
    gainSatFlag = false;
    
    fprintf('\n Try gain = ')
    for gind = 1:length(gainVals)
        fprintf('%d ', gainVals(gind))
        %-- Check that the previous gain setting wasn't saturated
        if gainSatFlag
            % If gainFlag is true, a prior gain setting was already
            % saturated so let's not bother measuring this gain setting
            
            % Let's provide a specific value in the data matrices though so
            % that we know this setting was considered
            FMTOmeas(pind, gind, :) = wierdVal2;
            PMmeas(pind, gind, :)   = wierdVal2;
            
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
            FMTOmeas(pind, gind, :) = wierdVal1;
            PMmeas(pind, gind, :)   = wierdVal1;
            
            fprintf('sat, ')
            
            % Set the saturation flag true so that we skip any higher gains
            gainSatFlag = true;
        else
            % The power was below 10 V and we can keep running as normal
            
            % Measure the redPM 
            readPM = VFN_PM100D_read(PM, pmNread);
            
            % Save the results
            FMTOmeas(pind, gind, :) = read;
            PMmeas(pind, gind, :)   = readPM;
        end
        
    end
    disp(' ')
end

%-- Remind user to set the laser power down again
power0reminder = 'DONE MEASURING::: You can now set the laser power low again';
fprintf('\n\n')
fprintf(power0reminder)
fprintf('\n\n')

%-- Set gain to lowest setting
FMTO.FMTO_scale = 5;
VFN_FMTO_setGain(FMTO);

%% cleanup all connections
VFN_cleanUpFMTO;
VFN_cleanUpPM100D;

%% Do some analysis
%-- Check if any matrix elements remained unassigned
if any(isnan(FMTOmeas(:)))
    error('Some element in FMTO data matrix is nan')
end
if any(isnan(PMmeas(:)))
    error('Some element in redPM data matrix is nan')
end

%-- Make a copy of the data for now
FMTOmeas_raw = FMTOmeas;
PMmeas_raw = PMmeas;

%-- Change any of the "wierdVal#' entries to nan
FMTOmeas(FMTOmeas == wierdVal1) = nan;
PMmeas(PMmeas == wierdVal1) = nan;
FMTOmeas(FMTOmeas == wierdVal2) = nan;
PMmeas(PMmeas == wierdVal2) = nan;

%-- Compute the STD along third dimension
FMTOmeas_std = std(FMTOmeas, [], 3);
PMmeas_std = std(PMmeas, [], 3);

%-- Take averages along third dimension 
FMTOmeas = mean(FMTOmeas, 3);
PMmeas = mean(PMmeas, 3);

%-- Display the final data
FMTOmeas
PMmeas
%% Save the data
if isSave
    fprintf('\n\nSAVING DATA TO: %s\n', svnm)
    save(svnm, 'gainVals', 'pwrVals', 'Delay', 'wierdVal1', 'wierdVal2', ...
        'FMTOmeas', 'PMmeas', 'FMTOmeas_raw', 'PMmeas_raw');
end
%% Remind user one last time to set the power low
warning(['REMINDER:', power0reminder])

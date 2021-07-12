%{
This script checks the femto voltage behavior when the diode saturates

It applies an increasing power ramp on the laser and measures the
corresponding voltage on the femto. The user provides a gain setting to use
for the femto. The script purposely saturates the diode and measures the
voltage beyond saturation to check how the reported voltage behaves.

This is because Dimitri pointed out that some detectors actually peak in
votlage at saturation and then beyond saturation, they start to report
decreasing voltages. If the femto does this, it could explain some of the
anomolous behavior I see in my donut scans.
%}

close all; clear all;

%% Define test conditions

% Define the gain value to test. This is the n in 10^n gain settings
gainVal = 8;     

% Define the laser power settings to try
pwrVals = 0.05:0.05:0.9;

% Define the time to wait after changing the power setting
Delay = 2.5;

% Provide filename for saving data
svnm = '/media/Data_Drive/VFN/TestbedData/210428_COV8/GainSatCal_Sample1';
% Flag marking whether to save data or not
isSave = false;

%% Femto setup
disp('Setting up femto')
% Set the number of reads to do in a single burst
FMTO.Nread   = 100;      % power samples in a single burst

% Connect to the femto
VFN_setUpFMTO;  

%% Laser setup
disp('Setting up laser')
% Connect to the laser
VFN_setUpKLS;

%% redPM setup
disp('Setting up redPM')
pmNread = 100;      % Number of samples to take
pmCalWvl = 635; % Wavelength for redPM for normalization

% Connect to the redPM
VFN_setUpPM100D

% Set operating wavelength for calibration
VFN_PM100D_setWavelength(PM, pmCalWvl);

% Wait a moment just in case the wavelength setting needs to settle
pause(1)

%% Take Data
%-- Preallocate data matrices
% Define matrix to hold femto measurements. 
   % Axes: (powerValue, samples in a single burst)
FMTOmeas = nan(length(pwrVals), FMTO.Nread);
% Define matrix to hold redPM measurements. 
   % Axes: (powerValue, samples in a single burst)
PMmeas = nan(length(pwrVals), pmNread);

%-- Iterate through power settings and take measurements
for pind = 1:length(pwrVals)
    % Set the laser power
    VFN_KLS_setPower(KLS,pwrVals(pind))
    % Wait for laser to stabilize, just to be safe
    pause(2)
    
    fprintf('\n Laser power set to: %f %%', pwrVals(pind))
    
    %-- Set the desired gain
    FMTO.FMTO_scale = gainVal;
    % Set isBlock=false so that we can use the user-provided delay here
    VFN_FMTO_setGain(FMTO,false);
    % Wait the desired time for the value to settle, just to be safe
    pause(Delay)
        
    %-- Take Measurements
    % Read power on the Femto
    read = VFN_FMTO_readN(FMTO);
        
    % Check that the power is less than 10
    if find(read>10)
        % The power was above 10. Warn the user
        fprintf('\nFemto is saturated: %f [V]\n',mean(read))
    end
    % Measure the redPM 
    readPM = VFN_PM100D_read(PM, pmNread);

    % Save the results
    FMTOmeas(pind, :) = read;
    PMmeas(pind, :)   = readPM;
end

%-- Set laser to lowest power
VFN_KLS_setPower(KLS,pwrVals(1))
%-- Set gain to lowest setting
FMTO.FMTO_scale = 5;
VFN_FMTO_setGain(FMTO);

%% cleanup all connections
VFN_cleanUpFMTO;
VFN_cleanUpKLS;
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

%-- Compute the STD along third dimension
FMTOmeas_std = std(FMTOmeas, [], 3);
PMmeas_std = std(PMmeas, [], 3);

%-- Take averages along third dimension 
FMTOmeas = mean(FMTOmeas, 2);
PMmeas = mean(PMmeas, 2);

%-- Display the final data
FMTOmeas
PMmeas

%% Plot the results
figure();
plot(PMmeas,FMTOmeas, '-o');
xlabel('redPM Power [W]')
ylabel('Femto Power [V]')
title(sprintf('Femto Saturation Response (gain = %d)',gainVal))
%% Save the data
if isSave
    fprintf('\n\nSAVING DATA TO: %s\n', svnm)
    save(svnm, 'gainVal', 'pwrVals', 'Delay', 'FMTOmeas', 'PMmeas', 'FMTOmeas_raw', 'PMmeas_raw');
end

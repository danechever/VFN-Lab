%{
Script for obtaining Femto-to-redPM Ratio - (using superK Ext-4)

*** NOTE: THIS SCRIPT ONLY WORKS WITH THE SUPERK EXTREME. I wasn't able to 
    easily get the comms working to the superK compact.

*** NOTE: The script only iterates through wavelengths for now. Could update
    it to do bandwidths but this is not needed at the moment.

This script uses the femto, redPM, and superK ext-4 to measure the femto
or redPM power at various wavelengths. These measurements can be used to
obtain the ratio between the femto voltage and redPM for a given wavelength.

The script procedure is:
- Select device on which to measure power (Femto or redPM)
- Set superK power
- Set superK bandwidth
- Set superK wavelength (lam)
- Measure power on selected device
  - If Femto, switch through all gain settings. 
- Set next superK lam. Iterate to sample all superK settings
- At the end, set the power and emission back to what they
  were at the start.

Optical setups should be:
- superK ext-4 output plugged into either the redPM or the femto

The output is a single matrix of power/voltage readings where each row 
represents a different lam setting and each column represents a gain 
setting. NOTE: since the redPM has no gain settings, redPM data matrices 
are just a column vector.

NOTE: the script avoids sampling gains that are/will be saturated.
%}


close all; clear all;
addpath(genpath(['..']));

%% Change to python2 environment (the NKT code only works from py2 for now)
pyVer = pyenv;
if ~strcmp(pyVer.Version, '2.7')
    % We aren't in the python2 (2.7) environment, let's fix that.
    path2py2env = '/home/vfndev/anaconda3/envs/py2_test/bin/python2';
    pyenv('Version', path2py2env);
    disp('Python Environment changed to python2')
end

%% Define test conditions

% Choose which measuremente device to use for power samples
measdev = 'redPM';   % must be either "Femto" or "redPM"

% Define number of sample points in a given burst of data
Nread   = 100;

% Define the gain values to try. Each value is the n in 10^n gain settings
  % Values in this vector should be in assending order
  % THIS ONLY APPLIES TO THE FEMTO
gainVals = [5, 6, 7, 8, 9, 10, 11];     

% Define the superK settings to use for the test
LASpower = 100;     % [%]  power to use (single value only)
LASLams  = [603,625,650,675,697]; % [nm] central wavelengths to sample
LASBW    = 6;       % [nm] Bandwidth to use (signle value only)

% Calibration wavelengths to use
  % Hypothetically this should be the same as "LASLams"
  % However, for the VFN poly paper, I forgot to change the redPM
    % calirbation wavelength during the wavelength-scan tests. As such, I
    % need to run this calibrator at a single cal wavelength on the redPM.
    % Setting this variable to a vector of contant values let's me do that.
redPM_calwvl = [603,603,603,603,603];

% Define the time to wait after changing the power/gain setting
Delay = 2.5;

% Define the value to use when a datapoint should be ignored
wierdVal1 = -8888;   % Value for when the current gain setting is saturated
wierdVal2 = -9999;   % Value for when a previous gain setting was saturated

% Provide filename for saving data
svnm = sprintf('/media/Data_Drive/VFN/TestbedData/210428_COV8/Femto2redPMCalibration_onSuperKEXT_LAM_%s_Sample13',measdev);
% Flag marking whether to save data or not
isSave = false;

%% Quick check 
%-- Doublecheck that the redPM calwvl is the same length as the lams
if length(redPM_calwvl) ~= length(LASLams)
    error('please provide equal length vectors for redPM calwvl and LASLams')
end
    
%% superK setup
disp('Setting up superK')
% Connect to the laser
VFN_setUpNKT;
% Make sure that the xtrmModuleAddr is set to the right one for the superK EXT
NKT.nktobj.xtrmModuleAddr = int64(15);

%-- Get starting power level
LASpower0 = VFN_NKT_getPowerLevel(NKT);

if LASpower0 ~= LASpower
    % Starting power is not what the user requested; set power
    fprintf('Starting power was: %0.2f%% ... Set to %0.2f%%\n',LASpower0,LASpower)
    VFN_NKT_setPowerLevel(NKT, LASpower);
else
    fprintf('Starting power was already %0.2f%%, did not change to %0.2f%%\n',LASpower0,LASpower)
end


%-- Get starting emission state
emsn0 = NKT.nktobj.get_emission();

%-- Set emission state 
if emsn0
    fprintf('SuperK was already on, not changing emission state')
else
    fprintf('superK emssion was off, setting emission on')
    VFN_NKT_setEmission(NKT, true);
end

%% Connect to the desired measurement device

if strcmpi(measdev, 'femto')
    %-- Femto setup
    disp('Setting up femto')
    % Set the number of reads to do in a single burst
    FMTO.Nread   = Nread;      % power samples in a single burst

    % Connect to the femto
    VFN_setUpFMTO;
elseif strcmpi(measdev, 'redpm')
    %-- redPM setup
    disp('Setting up redPM')

    % Connect to the redPM
    VFN_setUpPM100D

    % Set operating wavelength for calibration
    %VFN_PM100D_setWavelength(PM, LASLams(1));
    
    %-- The redPM doesn't need gain settings; set this vector to single nan
    gainVals = nan;
else
    error('The selected measurement device ("%s") was not recognized',measdev)
end

% Wait a moment just in case the wavelength setting needs to settle
pause(1)

%% Take Data
%-- Preallocate measurements data matrix
  % Axes: (Central wavelength setting, samples in a single burst, gainValues)
  % NOTE: 3rd dimension is implicitly omitted for redPM like this
meas = nan(length(LASLams), Nread, length(gainVals));

%-- Iterate through central wavelengths settings and take measurements
for lamind = 1:length(LASLams)
    %-- Set wavelengths on varia
    varlow = LASLams(lamind)-LASBW/2;
    varhi  = LASLams(lamind)+LASBW/2;
    VFN_NKT_setWvlRange(NKT, varlow, varhi);
    % Wait for superK to stabilize, just to be safe
    pause(2)
    
    fprintf('\n Varia set to: %0.2f x %0.2fBW [%0.2f - %0.2f]', LASLams(lamind), LASBW, varlow, varhi)

    if strcmpi(measdev, 'femto')
        %-- Deal with Femto measurement
        
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

                % Let's provide a specific value in the data matrix though 
                % so that we know this setting was considered
                meas(lamind, :, gind) = wierdVal2;

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
                meas(lamind, :, gind) = wierdVal1;

                fprintf('sat, ')

                % Set the saturation flag true so that we skip any higher gains
                gainSatFlag = true;
            else
                % The power was below 10 V and we can keep running as normal
                % Save the result
                meas(lamind, :, gind) = read;
            end
        end
    elseif strcmpi(measdev,'redpm')
        %-- Deal with redPM measurement
        
        % Set redPM operating wavelength for calibration
        VFN_PM100D_setWavelength(PM, redPM_calwvl(lamind));
        
        % Wait a moment before measurement just to be extra super safe
        pause(Delay)
        
        % Measure the redPM 
        read = VFN_PM100D_read(PM, Nread);
        
        % Save the result
        meas(lamind, :)   = read;        
    end
    disp(' ')
end

%-- Set superK to starting power and emission too
if LASpower0 ~= LASpower
    fprintf('superK set back to %0.2f%% power\n',LASpower0)
    VFN_NKT_setPowerLevel(NKT, LASpower0);
end
if ~emsn0
    fprintf('superK emssion set back to off')
    VFN_NKT_setEmission(NKT, emsn0);
end

%% cleanup all connections
%-- Disconnect from superK
VFN_cleanUpNKT;

%-- Disconnect from measurement device
if strcmpi(measdev,'femto')
    %-- Set gain to lowest setting
    FMTO.FMTO_scale = 5;
    VFN_FMTO_setGain(FMTO);
    %-- Close connection
    VFN_cleanUpFMTO;
elseif strcmpi(measdev,'redpm')
    VFN_cleanUpPM100D;
end

%% Do some analysis
%-- Check if any matrix elements remained unassigned
if any(isnan(meas(:)))
    error('Some element in data matrix is nan')
end

%-- Make a copy of the data for now
meas_raw = meas;

%-- Change any of the "wierdVal#' entries to nan
meas(meas == wierdVal1) = nan;
meas(meas == wierdVal2) = nan;

%-- Compute the STD along third dimension
meas_std = std(meas, [], 2);

%-- Take averages along third dimension 
meas = mean(meas, 2);

%-- Display the final data
meas
%% Save the data
if isSave
    fprintf('\n\nSAVING DATA TO: %s\n', svnm)
    save(svnm, 'gainVals', 'LASpower', 'LASBW', 'LASLams', 'redPM_calwvl', ...
        'Delay', 'wierdVal1', 'wierdVal2', 'meas', 'meas_raw', 'measdev');
end
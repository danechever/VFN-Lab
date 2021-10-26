%{ 
 Main_Transmissive
 Code for doing a normal scan with the VFN Transmissive bench
 - Uses the PIstages for fiber X, Y, and Z
 - Uses Zabers for vortex X and Y
 - Uses the redPM and a Zaber for in-beam normalization
 - Uses the Femto for through-fiber power measurements

** Based off VFN_Main_Transmissive

 * Allows normalization using the redPM or an MMF and the Femto
   - Also possible to not normalize at all
   * The MMF values reported in the fits files are -9999 if MMF is
        not used. This follows the same idea as the redPM when not used.
 * Outputs cubes of the same dimensionality as all previous code
   - ie: [X, Y, Nread, Focus, VortX, VortY]
 * Uses functions from VFN Analysis Library to interpret the data immediately
   - Like this, the values reported by this code are the same as we'll see when
    running the analysis code later (assuming same analysis params are used)
   - Full analysis (Fresnel losses, prop. losses, lens losses, etc.) is done.
 * Automatically prints values and results into the DataNotes file 
 * Uses the latest control libraries for the VFN devices
 * Calculates radially-averaged profile to determine proper planet coupling
 * Calculates rel. int. time using the null and the rad. avg. planet coupling 
   - Will optionally iterate through all frames and find the best rel. int.
 * Tries to not move the Zabers unless absolutely necessary.
 * Uses relative paths when possible (ex. for imports)

 Author: Daniel Echeverri (dechever@caltech.edu)
 Origin Date: 07/03/2020

%}


%% add the necessary functions
addpath(genpath(['..' filesep 'ControlCode']));
addpath(genpath(['..' filesep 'AnalysisCode' filesep 'AnalysisLib']));
addpath(genpath(['..' filesep '..' filesep 'VFN-Simulations' filesep 'VFNlib']));

close all
clear all

%% General settings 

% Directory for saving data:
svFld = '/media/Data_Drive/VFN/TestbedData/211019_NLS1';

%-- Experiment name for figures
expNm = '47Don_635VortScan13';

%-- Custom auto save text
% Change at every run
run_msg = 'Set foc back to that from #45 then scan foc and vortex in X and Y.';
% Change when doing new set of experiments
set_msg = 'Optimize vortex centering.';
las_msg = 'KLS635, 5.0mW.';   % Laser/Varia settings
flt_msg = 'N/A';                     % Optical Filter used
vtx_msg = 'JPL Poly V550-BB15 Charge 2';            % Vortex in use

%-- Final analysis stuff
isRadAvgAnalysis = true;    % calc rel. int. time in all frames using rad avg

%~~ ZABER STUFF 
% Flag to enable vortex scan.
isVortScan = true;

% Distance to move zabers for backlash removal [in mm]
zabBacklash = 0.040;         %Note: affects vortex

% Vortex center of scan [in mm]
VXcenter = 15.5068;       % 
VYcenter = 12.8504;       %13.828012-6 (5.8410 is min if X is at 2.8444)

% Vortex scan properties
%To include center values as a point in scan; use ODD number of points
VXpoints = 3;           % Number of X points in scan
VYpoints = VXpoints;    % Number of Y points in scan

% Vortex step params in [mm]  
vStepRange = 0.04; %Vortex will be scanned +/- this value
                    % ie. scan will be: 
                    % [VXcen-vStepRa to VXcen+vStepRa] w/ Vxpoi steps
%~~ END ZABER STUFF 

%~~ PI STUFF 
% Fiber Focus Scan properties
%To include center value as a point; use ODD number of points
Zcenter   = 3.400000;% MMF: 9.193358;  % [mm] Focus (Z) center point
Zpoints   = 3;           % Number of focci taken (exact number; no longer +1)
ZStepSize = 0.010;      % Step size in mm

% Fiber X/Y center of scan in mm
Xcenter = 8.477166;% MMF: 7.405309;  % [mm]
Ycenter = 8.208833;% MMF: -9.856926; % [mm]

% Fiber scan properties
%To include center values as a point in scan; use EVEN number of points
Xpoints = 30;% number of X points in scan (actual scan will have +1)
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle

% Fiber step sizes in Volts
refStep   = 10; % refStep is the step size for a full PSF with a 10*10 grid
StepSize  = refStep/(Xpoints/1e-3);
%~~ END PIEZO STUFF 

%~~ FEMTO POWER METER STUFF 
FMTO.isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO.FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
FMTO.Nread   = 100;      % power samples at a given locaiton
%Nrate   = 1000;     % rate of scan [samples/second]
% Add delay for settling time of detector
% NOTE::: Chose arbitrary delay values for now
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)
%~~ END FEMTO POWER METER STUFF 

%~~ FEMTO NORMALIZATION STUFF
isMMFNorm = false;
%~~ END FEMTO NORMALIZATION

%~~ RED (NORM) POWER METER STUFF 
pmNread = 100;      % Number of samples to take
isPMNorm = true;   % Flag to mark whether Thorlabs PM100D should be read
                        % This flag is useful in case the PM is not in the
                        % system. -9999 will replace the pmRead values.
% pmNormReport replaced with nrmValReport to reflect the fact that it is
    % representing the nrmVal variable, not pmRead1 anymore.
nrmValReport = 1;   % @635nm typ = 14.08; @780nm typ = 44.1
                        % Valu for fibertip power when isPMNorm=false. Only
                        % affects the printed values at end; NORM'd cube
                        % values are still un-normed. SET =1 to not norm
                        % the printed values.
pmCalWvl = 635; % Wavelength for redPM for normalization
%~~ END RED (NORM) POWER METER STUFF

%% Zaber Setup
fprintf('---Performing Zaber Setup\n')
isZab = false;
if isVortScan || isPMNorm
    % Enable zabers since they will be used at some point
    isZab = true;
end

%~~ Connect to Zabers to query position
VFN_setUpBenchZabers; % Instantiate and rename the zabers

%- Make sure pmX is out of the way
if isfield(Zabs, 'pmX')
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
elseif isPMNorm
    error('PM Norm is requested but no pmX Zaber exists.')
end

%- Query positions
if (VFN_Zab_getPos(Zabs.vortX) ~= VXcenter) || (VFN_Zab_getPos(Zabs.vortY) ~= VYcenter)
    % vort not at cent. pos; if not doing vort scan, enable zabs to correct pos later
    isZab = true;
end

%- Disable Zabers for rest of code if not needed anymore
if ~isZab 
    VFN_cleanUpZabers;
end
    

%% PI stage (nulling fiber stage) Setup
fprintf('---Performing PI Stage Setup\n')
% Connect to the MDT Piezos
VFN_setUpBenchPIStages;

%% Femto setup
fprintf('---Performing Femto Setup\n')
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
FMTO.gnFact  = 10.01;

% Setup Femto 
VFN_setUpFMTO;  

% Set FMTO gain to user-provided value (known start value for code)
VFN_FMTO_setGain(FMTO);

fprintf('Current Gain setting: %i\n', FMTO.FMTO_scale)

%% Red Thorlabs PM setup
fprintf('---Performing redPM Setup (if needed)\n')
%-- Connect to red PM only if it will be used
if isPMNorm
    VFN_setUpPM100D
    
    % Set operating wavelength for calibration
    VFN_PM100D_setWavelength(PM, pmCalWvl);
end

%% DataNotes automation pt1
% Write the filename
datName = [filesep 'DataNotes.txt'];
datFl = strcat(svFld, datName);
fileID = fopen(datFl, 'a');
fprintf(fileID, '\r\n- %s\r\n', expNm);

% General comment
fprintf(fileID, '    %s\r\n', run_msg);
fprintf(fileID, '    * %s\r\n', set_msg);

% Scan Params section
fprintf(fileID, '    Scan params:\r\n');
fprintf(fileID, '        isZab:         %s\r\n', mat2str(isZab));
fprintf(fileID, '        isVortScan:    %s\r\n', mat2str(isVortScan));%    
trash = [VXcenter; VYcenter];
fprintf(fileID, '        Vortex pos:    x = %9.6fmm    y = %9.6fmm\r\n', trash);
fprintf(fileID, '        Vbacklash:     %5.3fmm\r\n', zabBacklash); 
trash = [VYpoints; VXpoints; vStepRange];
fprintf(fileID, '        v scan param:  vYpoi = %d, vXpoi = %d, vStepRange = %4.2f\r\n', trash);
fprintf(fileID, '        Laser power:   %s\r\n', las_msg);
fprintf(fileID, '        Optical Filt:  %s\r\n', flt_msg);
fprintf(fileID, '        Vortex Mask:   %s\r\n', vtx_msg);
trash = [Xcenter; Ycenter; Zcenter];
fprintf(fileID, '        scan center    (%8.6f mm,%8.6f mm, %8.6f mm)\r\n', trash);
trash = [Xpoints; Zpoints; refStep; refStep; Xpoints; ZStepSize];
fprintf(fileID, '        Scan params:   xpoi = %d, zpoi = %d, refStep = %d, stepsize = %d/(%d/1e-3); zstep = %5.3f mm\r\n', trash);
fprintf(fileID, '        isAutoScale:   %s\r\n', mat2str(FMTO.isAutoScale));
trash = [FMTO.FMTO_scale; FMTO.Nread];
fprintf(fileID, '        Femto:         FMTO_scale = %d; Nread = %d\r\n', trash);
trash = [nrmValReport; pmCalWvl];
fprintf(fileID, '        redPM norm:    isPMNorm = %s; nrmValReport = %5.2f; pmCalWavelength = %dnm\r\n', mat2str(isPMNorm), trash);
fprintf(fileID, '        isMMFNorm: %s\r\n\r\n', mat2str(isMMFNorm));
fclose(fileID);

%~~ END DataNotes pt1 automation

%% Prepare scan
%-- Define Fiber scan locations
% Colon indexing for fib X,Y so that even xPoi has xCent as middle point in scan
distX=(Xcenter-StepSize*Xpoints/2):StepSize:(Xcenter+StepSize*Xpoints/2);
distY=(Ycenter-StepSize*Ypoints/2):StepSize:(Ycenter+StepSize*Ypoints/2);
% Linspace for foc so that number of points is exactly what user requested, not +1
distZ = linspace(Zcenter-ZStepSize*floor(Zpoints/2),Zcenter+ZStepSize*floor(Zpoints/2),Zpoints);

%-- Define Vortex scan locations
% Pre-allocate distVX and distVY for case when they shouldn't be changed
distVX = VXcenter;
distVY = VYcenter;

if isVortScan && (VXpoints > 1)
    % Change distVX since it is being scanned
    distVX = linspace(VXcenter-vStepRange, VXcenter+vStepRange, VXpoints);
end
if isVortScan && VYpoints > 1
    % Change distVY since it is being scanned
    distVY = linspace(VYcenter-vStepRange, VYcenter+vStepRange, VYpoints);
end

%-- Preallocate data matrices
meas= nan(length(distX),length(distY),FMTO.Nread,length(distZ),length(distVX), length(distVY)); 
measScl = meas;         %Matrix for semi-processed data
sclSz   = size(meas); sclSz(3) = 3;
scales  = nan(sclSz);   %Matrix for scales, gains, and biases
%fibPos = scales;        %Matrix for fiber positions (x,y,z)

%% Check that stage limits will not be exceeded
%-- Check zaber limits
if (distVX(1)-zabBacklash < Zabs.vortX_lower) || (distVX(end) > Zabs.vortX_upper)
    error('vortX motion will exceed limits')
end
if (distVY(1)-zabBacklash < Zabs.vortY_lower) || (distVY(end) > Zabs.vortY_upper)
    error('vortY motion will exceed limits')
end


%-- Check PI limits
if (distX(1) < PIdevs.fibX_lower) || (distX(end) > PIdevs.fibX_upper)
    error('fibX motion will exceed limits')
end
if (distY(1) < PIdevs.fibY_lower) || (distY(end) > PIdevs.fibY_upper)
    error('fibY motion will exceed limits')
end
if (distZ(1) < PIdevs.fibZ_lower) || (distZ(end) > PIdevs.fibZ_upper)
    error('fibZ motion will exceed limits')
end


%% Perform scan

if isZab && (distVX(1) ~= VFN_Zab_getPos(Zabs.vortX))
	%remove zaber (X) backlash
    VFN_Zab_move(Zabs.vortX, distVX(1)-zabBacklash);
end
for a = 1:length(distVX)
    if isZab && (distVX(a) ~= VFN_Zab_getPos(Zabs.vortX))
        %move vortex to new X position
        fprintf('\nVortX pos: %f',VFN_Zab_move(Zabs.vortX, distVX(a)));
    end

    %remove zaber (Y) backlash
    if isZab && (distVY(1) ~= VFN_Zab_getPos(Zabs.vortY))
        % Remove y backlash now that motion is confirmed safe
        VFN_Zab_move(Zabs.vortY, distVY(1)-zabBacklash);
    end
    
    for b = 1:length(distVY)
        if isZab && (distVY(b) ~= VFN_Zab_getPos(Zabs.vortY))
            %move vortex to new Y position
            fprintf('  VortY pos: %f\n',VFN_Zab_move(Zabs.vortY, distVY(b)));
        end

        for k=1:length(distZ)
            % Move fibZ
            VFN_PIStage_move(PIdevs.fibZ, distZ(k));
                        
            for j=1:Xpoints+1
                % Print current scan status
                fprintf('Vort %i of %i | Foc %i of %i | Column %i of %i\n', ...
                        (a-1)*length(distVY)+b,length(distVY)*length(distVX),...
                        k, length(distZ), j, length(distX))
                
                % Move fibX
                VFN_PIStage_move(PIdevs.fibX, distX(j));
                                
                for i=1:Ypoints+1
                    % Move fibY
                    VFN_PIStage_move(PIdevs.fibY, distY(i));
                    
                    % Let femto voltage settle at new position
                    pause(Delay);
                                       
                    % Take a sample reading at current power
                    read    = VFN_FMTO_readN(FMTO);
                    
                    % Perform auto scaling if desired
                    if FMTO.isAutoScale
                        % Save old FMTO_scale for comparison after autoGain
                        old_scale = FMTO.FMTO_scale;
                        % Modify gain accordingly
                        FMTO = VFN_FMTO_setAutoGain(FMTO, mean(read));
                        % Check if re-read is needed
                        if old_scale ~= FMTO.FMTO_scale
                            % FMTO_scale changed so the gain changed
                            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO.FMTO_scale)
                            read    = VFN_FMTO_readN(FMTO);
                            if FMTO.FMTO_scale > 9
                                warning('Gain >9: %i',FMTO.FMTO_scale)
                            end
                        end
                    end
                    if find(read>10)
                        warning('Power is too high')
                    end
                    
                    % Save fiber positions
                    %fibPos(j,i,1) = VFN_PIStage_getPos(PIdevs.fibX);
                    %fibPos(j,i,2) = VFN_PIStage_getPos(PIdevs.fibY);
                    %fibPos(j,i,3) = VFN_PIStage_getPos(PIdevs.fibZ);

                    % Store data
                    meas(j,i,:,k,a,b) = read;
                    scales(j,i,1,k,a,b) = FMTO.FMTO_scale;
                    scales(j,i,2,k,a,b) = FMTO.gnFact^-(FMTO.FMTO_scale-6);
                    
                    % Get Bias at current gain setting
                    locBias = VFN_FMTO_getBias(FMTO.FMTO_scale);
                    
                    % Store bias data
                    scales(j,i,3,k,a,b) = locBias;
                    measScl(j,i,:,k,a,b) = read - locBias; %Subtract bias for semi-processed matrix
                end
            end
        end
    end
end

%% Measure power at Fiber Tip
%-- Use the MMF to measure the power
    % Do MMF first so that redPM is gauranteed to be out of beam for measurement
if isMMFNorm
    disp('Switch fibers. Press any key to continue')
    %%%% When this becomes automated, this will become a PI command
    %%%% rather than a manual switch CHANGE_MARK
    pause;
    %Now the MMF should be in the beam and can collect all the light
    %Use the femto to collect the total power, same as collecting data
    mmf_read    = VFN_FMTO_readN(FMTO);
    % Perform auto scaling if desired
    if FMTO.isAutoScale
        % Save old FMTO_scale for comparison after autoGain
        old_scale = FMTO.FMTO_scale;
        % Modify gain accordingly
        FMTO = VFN_FMTO_setAutoGain(FMTO, mean(mmf_read));
        % Check if re-read is needed
        if old_scale ~= FMTO.FMTO_scale
            % FMTO_scale changed so the gain changed
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO.FMTO_scale)
            mmf_read    = VFN_FMTO_readN(FMTO);
            if FMTO.FMTO_scale > 9
                warning('Gain >9: %i',FMTO.FMTO_scale)
            end
        end
    end
    if find(mmf_read>10)
        warning('Power is too high')
    end
    
    % Record femto gain setting during MMF read
    MMF_FMTO_scale = FMTO.FMTO_scale;    
else
    %-- Set mmf_read and MMF_FMTO_scale to -9999 if MMF is not used
    mmf_read = ones(FMTO.Nread,1)*-9999;
    MMF_FMTO_scale = -9999;
end

%-- Use the red power meter to measure the power
if isPMNorm
    %-- Move redPM into the beam
    fprintf('\nredPM norm used. Moving fibZ...')
    
    %Move fiber stage out of the way to avoid collision with pmX
    VFN_PIStage_move(PIdevs.fibZ, PIdevs.fibZ_lower);
    fprintf('  fibZ at %0.2f',VFN_PIStage_getPos(PIdevs.fibZ))
    
    % Brief pause just to be safe
    pause(2);
    
    %Move if no collision
    fprintf('\nMoving pmX...')
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_upper);
    fprintf(' Doing measurments...')

    %-- Read from red PM 
    % Wait for PM value to settle
    pause(2);

    % Take measurements of the power (use try catch to avoid errors)
    NTries = 5;
    for i = 1:NTries
        try 
            % This is where the read often fails 
            pmRead = VFN_PM100D_read(PM, pmNread);
        catch ME
            % An error occurred, let the user know
            warning('PM read failed')
            if i < 5
                % We have tried less than NTries times, try again
                
                % Wait a moment to give the redPM time to figure itself out
                pause(1);
                
                % First, try playing with the wavelength settings
                    % This seems to help for some reason
                VFN_PM100D_setWavelength(PM, pmCalWvl);
                pause(0.1);
                VFN_PM100D_getWavelength(PM);
                
                % Now try again (continue will jump straight to next iter.)
                fprintf('Trying again\n')
                continue
            else
                % We've tried 5 times, it's time to give up
                error('Multiple redPM reads failed')
            end
        end
        % If we get here, succeeded in reading so exit for-loop and proceed
        break
    end
    
    % Move redPM zaber, pmX, out of beam again
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
    fprintf('\nDone. pmX moved to %0.2f',VFN_Zab_getPos(Zabs.pmX))
else
    %-- Set pmRead to -9999 if PM is not to be read
    pmRead = ones(pmNread,1)*-9999;
end

%% close the connections to all devices
% Clean up FMTO
VFN_cleanUpFMTO;
fprintf('\n Femto Gain set to %i',FMTO.FMTO_scale);

if isZab
    %-- Return vortex to center position if needed
    if VXcenter ~= VFN_Zab_getPos(Zabs.vortX)
        % Remove Zab (X) backlash
        VFN_Zab_move(Zabs.vortX, VXcenter-zabBacklash);
        % Center Zab (X)
        VFN_Zab_move(Zabs.vortX, VXcenter);
    end
    if VYcenter ~= VFN_Zab_getPos(Zabs.vortY)
        % Remove Zab (Y) backlash
        VFN_Zab_move(Zabs.vortY, VYcenter-zabBacklash);
        % Center Zab (Y)
    	VFN_Zab_move(Zabs.vortY, VYcenter);
    end
    
    % Print position of vortex
    fprintf('\nVortX pos: %f  VortY pos: %f',...
        VFN_Zab_getPos(Zabs.vortX), ...
        VFN_Zab_getPos(Zabs.vortY));
    
    
    %-- Move Calibration PM out of the beam (again; just to be safe)
    if isPMNorm
        fprintf('  |  PM Zab pos: %f\n', VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower));
    else
        fprintf('\n');  % terminate line from zaber recentering
    end
    
    %-- Clean up
    VFN_cleanUpZabers;
end

%Set PI stages to central position
VFN_PIStage_move(PIdevs.fibX, Xcenter);
VFN_PIStage_move(PIdevs.fibY, Ycenter);
fprintf('\nfibZ pos: %f', VFN_PIStage_move(PIdevs.fibZ, Zcenter));
% Close connection to the PI Stages
VFN_cleanUpPIStages

%-- Disconnect from calibration PM
if isPMNorm
    %-- Only disconnect if we connected earlier
    VFN_cleanUpPM100D;
end

%% Scale data accordingly to compensate for varied gain
%-- This is for plotting and analysis; the raw data is still stored as well

%** Bias is already subtracted in the main for loop

%-- Define the propagation and fresnel reflection loss values
FL = 0.9654;  % (1-loss) Fraction per surface
PL = 0.9966;  % (1-loss) Fraction per meter of fiber

%-- Define scaled matrix; equal to mean of meas unless gain was varied
measScl = mean(measScl,3).*scales(:,:,2,:,:,:);

%-- Get FMTO sensitivity at given wavelength
%fmtoSnsScl = VFN_getFmtoResponsivity(pmCalWvl);
fmtoSnsScl = 0.62; % value from manual comparison

%-- Define a separate matrix which normalizes the data by power on the 
if isPMNorm
    %-- When actual PM reading was done, normalize as usual
    % Account for loss from fiber lens
        % 3.4% reflection spec at 780nm, 0.3% reflection measured at 635nm
    %pmRead1 = mean(pmRead)*(1-0.034);%0.997;    % Updated on 2/8/19 based on 10/26/18
    % Convert the calibration PM value from W to uW
    pmRead1 = mean(pmRead)*10^6;
    % Translate the calibration PM value from uW to V (at 10^6)
    pmRead1 = pmRead1*fmtoSnsScl;     % V/uW conversion is ~ 0.63 @633; ~0.9 @ 780
    % Account for Fresnel and Propagation losses
    nrmVal = (1/pmRead1)*(1/(FL^2))*(1/(PL^2));
else
    %-- When PM was not read, set normalization value to 1
    nrmVal = 1;
end

if isMMFNorm
    % Get Bias at current gain setting
    locBias = VFN_FMTO_getBias(MMF_FMTO_scale);
    mmf_norm = mmf_read-locBias;
    % Calculate gain factor and rescale accordingly
    scl2 = FMTO.gnFact^-(MMF_FMTO_scale-6);
    mmf_norm = mean(mmf_norm).*scl2;
    
    % Account only for propagation losses since fresnel losses will be present
        % in both (SMF and MMF) fibers hence implicitly accounted for
    nrmVal = (1/mmf_norm)*(1/(PL^2));
end

% Create the matrix and normalize by total power on fib tip
measScl2 = measScl*nrmVal;


%% Display of data
for a = 1:length(distVX)
    for b = 1:length(distVY)
        for k = 1:length(distZ)
            figure();
            colormap('gray');
            %Not sure why transpose and "axis xy" are needed for correct orient.
            imagesc((distX-Xcenter)*1e3,(distY-Ycenter)*1e3,transpose(measScl2(:,:,:,k,a,b)));
            titStr = sprintf('PWR Thru Fib, Z=%0.5f mm, vort=(%0.5f, %0.5f)', ...
                                distZ(k), distVX(a), distVY(b));
            title(titStr)
            axis xy;axis image%axis square;
            cbr = colorbar;
            ylabel(cbr, 'Normed Coupling')
            set(gca,'TickDir','out')
            xlabel('Fiber X position [um]')
            ylabel('Fiber Y position [um]')
            tag = sprintf('%ix%i',Xpoints,Ypoints);
            if k ~= 1
                tag = sprintf([tag '_Foc%i'],k);
            end
            if a ~= 1 || b ~= 1
                tag = [tag sprintf('_Vort%ix%i',a,b)]; %#ok<AGROW>
            end
            saveas(gcf, [svFld filesep expNm '_' tag '_CoupMap.png'])
        end
    end
end

%% fits save
import matlab.io.*
tag = sprintf('%ix%i',Xpoints,Ypoints);
if Zpoints ~= 1
    tag = [tag sprintf('x%i_Cube',Zpoints)];
end
if isVortScan
    tag = [tag sprintf('_%ix%iVortScan', length(distVX), length(distVY))];
end

% Cube names
datCubeNm   = [svFld filesep expNm '_' tag '.fits'];
GnCubeNm    = [svFld filesep expNm '_' tag '_GAINS.fits'];
SRCubeNm    = [svFld filesep expNm '_' tag '_SEMIREDU.fits'];
NormCubeNm  = [svFld filesep expNm '_' tag '_NORMREDU.fits'];

% Define permute order:
permOrd = [2,1,3,4,5,6];

%-- SAVE DATA CUBE
fitmap= fits.createFile(datCubeNm);
fits.createImg(fitmap,'double',size(permute(meas,permOrd)));
%header data
fits.writeKey(fitmap, 'Xcenter ', Xcenter);
fits.writeKey(fitmap, 'Xpoints ', Xpoints);
fits.writeKey(fitmap, 'Ycenter ', Ycenter);
fits.writeKey(fitmap, 'Ypoints ', Ypoints);
fits.writeKey(fitmap, 'Zcenter ', Zcenter);
%fits.writeKey(fitmap,'Zstart ',Zstart);
fits.writeKey(fitmap, 'Zpoints ', Zpoints);
fits.writeKey(fitmap, 'XYsteps ', StepSize);
fits.writeKey(fitmap, 'Zsteps ',  ZStepSize);
fits.writeKey(fitmap, 'Nread ',  FMTO.Nread);
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', 'Nread');
fits.writeKey(fitmap, 'NAX4', 'Focus');
fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');
%Zaber status
fits.writeKey(fitmap, 'isZab ',  logical(isZab));
fits.writeKey(fitmap, 'VXcenter',  VXcenter);
fits.writeKey(fitmap, 'VYcenter',  VYcenter);
fits.writeKey(fitmap, 'isVortSc', logical(isVortScan));
fits.writeKey(fitmap, 'VXpoints', VXpoints);
fits.writeKey(fitmap, 'VYpoints', VYpoints);
fits.writeKey(fitmap, 'vStepRng', vStepRange);
% Power meter level
fits.writeKey(fitmap, 'isAutScl', logical(FMTO.isAutoScale));
ind = strfind(GnCubeNm ,filesep);
fits.writeKey(fitmap, 'DatCube ', GnCubeNm(ind(end-1)+1:end));
ind = strfind(SRCubeNm ,filesep);
fits.writeKey(fitmap, 'SRCube ', SRCubeNm(ind(end-1)+1:end));
ind = strfind(NormCubeNm ,filesep);
fits.writeKey(fitmap, 'NormCube ', NormCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NormFlag', logical(isPMNorm));
fits.writeKey(fitmap, 'MMFFlag', logical(isMMFNorm));
fits.writeKey(fitmap, 'NormMean', mean(pmRead));
fits.writeKey(fitmap, 'NormSTD', std(pmRead));
fits.writeKey(fitmap, 'MMFMean', mean(mmf_read));
fits.writeKey(fitmap, 'MMFSTD', std(mmf_read));
fits.writeKey(fitmap, 'MMFGain', MMF_FMTO_scale);
fits.writeKey(fitmap, 'FresLoss', FL);
fits.writeKey(fitmap, 'PropLoss', PL);
fits.writeKey(fitmap, 'NormWvl', pmCalWvl);
fits.writeKey(fitmap, 'NormScl', fmtoSnsScl);
fits.writeKey(fitmap, 'CNTRLCD', 'Main_Trans_PIStg');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(meas,permOrd));
fits.closeFile(fitmap);

%-- SAVE GAIN SETTINGS
fitmap= fits.createFile(GnCubeNm);
fits.createImg(fitmap,'double',size(permute(scales, permOrd)));
%header data
ind = strfind(datCubeNm,filesep);
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', '1=FMTO_Scale;2=Gain;3=Bias');
fits.writeKey(fitmap, 'NAX4', 'Focus');
fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(scales,permOrd));
fits.closeFile(fitmap);

%-- SAVE SEMI REDUCED DATA
fitmap= fits.createFile(SRCubeNm);
fits.createImg(fitmap,'double',size(permute(measScl, permOrd)));
%header data
ind = strfind(datCubeNm,filesep);
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', 'mean(meas-bias)*Gain');
fits.writeKey(fitmap, 'NAX4', 'Focus');
fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(measScl,permOrd));
fits.closeFile(fitmap);

%-- SAVE SEMI REDUCED NORM'D DATA
fitmap= fits.createFile(NormCubeNm);
fits.createImg(fitmap,'double',size(permute(measScl2,permOrd)));
%header data
ind = strfind(datCubeNm,filesep);
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', 'mean(meas-bias)*Gain/pmRead/FL/PL');
fits.writeKey(fitmap, 'NAX4', 'Focus');
fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');
fits.writeKey(fitmap, 'NormFlag', logical(isPMNorm));
fits.writeKey(fitmap, 'MFFlag', logical(isMMFNorm));
fits.writeKey(fitmap, 'NormMean', mean(pmRead));
fits.writeKey(fitmap, 'NormSTD', std(pmRead));
fits.writeKey(fitmap, 'MMFMean', mean(mmf_read));
fits.writeKey(fitmap, 'MMFSTD', std(mmf_read));
fits.writeKey(fitmap, 'MMFGain', MMF_FMTO_scale);
fits.writeKey(fitmap, 'FresLoss', FL);
fits.writeKey(fitmap, 'PropLoss', PL);

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(measScl2,permOrd));
fits.closeFile(fitmap);

%% Report results (basic analysis)

%-- Use user-provided norm value if desired
if ~(isPMNorm || isMMFNorm)
    nrmVal = nrmValReport;
end

%-- Find null in central region
cropVal = 7;
[mnVal, mnInd] = VFN_An_getEta_s(measScl, cropVal);
fprintf('\nMin in center region (un-normed): %f',mnVal)
fprintf(['\n Min Location: ', ...
         '\n  X index = %3i, X   =  %09.6f mm', ...
         '\n  Y index = %3i, Y   =  %09.6f mm', ...
         '\n  Foc ind = %3i, Foc =  %09.6f mm', ...
         '\n  VX ind  = %3i, VX  =  %09.6f mm', ...
         '\n  VY ind  = %3i, VY  =  %09.6f mm'], ...
         mnInd(1), distX(mnInd(1)), mnInd(2), distY(mnInd(2)), ...
         mnInd(4), distZ(mnInd(4)), mnInd(5), distVX(mnInd(5)), ...
         mnInd(6), distVY(mnInd(6)));


%-- Find overall max
[mxVal2, mxInd2] = VFN_An_getEta_p(measScl);
fprintf('\nOverall Max (un-normed):          %f',mxVal2)
fprintf('\n Overall Max Location: (%i, %i, %i, %i, %i, %i)', ...
        mxInd2(1), mxInd2(2), mxInd2(3), mxInd2(4), mxInd2(5), mxInd2(6));
    
%-- Find max in frame with min value
meas3 = measScl(:,:, mnInd(3), mnInd(4), mnInd(5), mnInd(6)); %Frame with min
mxVal = max(meas3(:)); 
fprintf('\nMax in min Frame (un-normed):     %f', mxVal)
fprintf('\nMax/Min, min Frame (un-normed):   %f', mxVal/mnVal)
if ~(isPMNorm || isMMFNorm)
    fprintf('\nWARN: redPM not used; fib tip power provided by user');
end

%-- Get normed values (used measScl to get only gain-corrected values above)
eta_s = mnVal*nrmVal;
eta_p = mxVal*nrmVal;

fprintf('\nPower on Fiber tip:      %f', 1/nrmVal);
fprintf('\nMin, normed (eta_s):     %f', eta_s);
fprintf('\nMax, normed (eta_p):     %f\n', eta_p);

if ~sum(strfind(expNm, 'PSF'))
    % Assume we are analyzing a donut
    
    %-- Find eta_p as radially-averaged max (in min frame)
        % meas3*nrmVal to get fully normed value directly
    [radProf, radVec] = VFN_An_radAverage(meas3*nrmVal,[mnInd(2), mnInd(1)]);
    radMax = max(radProf);
    fprintf('Max, normed (radial average): %f\n', radMax);

    %-- Display radially-averaged coupling line profile
    figure();
    z = plot(radVec,radProf*100, 'LineWidth', 3);
    title('Average Radial Coupling (min Frame)')
    ylabel('Coupling [%]');

    %-- Print pertinent values in command window
    fprintf('\n- Max/Min ratio in Frame = (%f/%f)    ------> = %0.2f\n', mxVal, mnVal, eta_p/eta_s)
    fprintf('- Eta_s = (%f*%f) w/ approx vals for thrpt -> = %f\n', mnVal, nrmVal, eta_s)
    fprintf('- Eta_p = (%f*%f) w/ approx vals for thrpt -> = %f\n', mxVal, nrmVal, eta_p)
    fprintf('- Average Eta_p w/ approx vals for thrpt          ------> = %f\n', radMax) 
	fprintf('- Rel. Int. Time = (%f/%f ^2) w/ rad Avg ---> = %f\n', eta_s, radMax, eta_s/(radMax^2))
    
    %-- Automated DataNotes
    fileID = fopen(datFl, 'a');
    fprintf(fileID, '\n      - Max/Min ratio in Frame = (%f/%f)    ------> = %0.2f\n', mxVal, mnVal, eta_p/eta_s);
    fprintf(fileID, '      - Eta_s = (%f*%f) w/ approx vals for thrpt -> = %f\n', mnVal, nrmVal, eta_s);
    fprintf(fileID, '      - Eta_p = (%f*%f) w/ approx vals for thrpt -> = %f\n', mxVal, nrmVal, eta_p);
    fprintf(fileID, '      - Average Eta_p w/ approx vals for thrpt          ------> = %f\n', radMax);
	fprintf(fileID, '      - Rel. Int. Time = (%f/%f ^2) w/ rad Avg ---> = %f\n', eta_s, radMax, eta_s/(radMax^2));
    fclose(fileID);
else
    % Analysing the PSF
    fprintf(['\n Max Location: ', ...
        '\n  X index = %3i, X   =   %0.6f     mm', ...
        '\n  Y index = %3i, Y   =   %0.6f     mm', ...
        '\n  Foc ind = %3i, Foc =   %0.6f  mm', ...
        '\n  VX ind  = %3i, VX  =   %0.6f mm', ...
        '\n  VY ind  = %3i, VY  =   %0.6f mm'], ...
        mxInd2(1), distX(mxInd2(1)), mxInd2(2), distY(mxInd2(2)), ...
        mxInd2(4), distZ(mxInd2(4)), mxInd2(5), distVX(mxInd2(5)), mxInd2(6), distVY(mxInd2(6)));
    
    fprintf('\n\n- Max (un-normed):            %f', mxVal2)
    fprintf('\n- Power on Fiber tip:         %f', 1/nrmVal)
    fprintf('\n- Max, normed (coupling):     %f\n', mxVal2*nrmVal)
    % Automated DataNotes
    fileID = fopen(datFl, 'a');
    fprintf(fileID, '\r\n         - Max:                        %f', mxVal2);
    fprintf(fileID, '\r\n         - Power on Fiber tip:         %f', 1/nrmVal);
    fprintf(fileID, '\r\n         - Max, normed (coupling):     %f\r\n\r\n', mxVal2*nrmVal);
    fclose(fileID);
    % Display max in each frame of focus
    tmp1 = nan(size(measScl2,4),1);
    for iii = 1:size(measScl2, 4)
        tmp1(iii) = max(max(measScl2(:,:,:,iii)));
    end
    fprintf('\n Coupling at all focci:\n')
    disp(tmp1')
end

%% Report Results (radial average analysis)
if isRadAvgAnalysis
%-- Preallocate data matrices 
% radially-averaged maximum in each frame
radMxs = nan(length(distZ), length(distVX), length(distVY));
% Null value in each frame
radMns = radMxs;

%-- Iterate through frames and calculate ratio
for a = 1:length(distVX)
    for b = 1:length(distVY)
        for k = 1:length(distZ)
            % Extract desired fram
            meas2 = measScl2(:,:,:, k, a, b);
            % Find null in central region
            [radMns(k, a, b), radMnInd] = VFN_An_getEta_s(meas2, cropVal);
            % Calculate average radial profile
            [radProf2, ~] = VFN_An_radAverage(meas2,[radMnInd(2), radMnInd(1)]);
            radMxs(k, a, b) = max(radProf2);                
        end
    end
end

%-- Calculate relative integration time in each frame
relTints = radMns./(radMxs.^2);

%-- Find best value
[bestTint, bestTintInd] = min(relTints(:)); 
[bestZ, bestVX, bestVY] = ind2sub(size(relTints), bestTintInd);

%-- Print results
fprintf('\nBest Integration Time: %f\n', bestTint);
fprintf('Best Location: (fibZ %i, vortX %i, vortY %i)\n', bestZ, bestVX, bestVY);
% Print into DataNotes
fileID = fopen(datFl, 'a');
fprintf(fileID, '      - Best Rel. Int. Time = (%f/%f ^2) w/ rad avg -> = %f\n',...
    radMns(bestZ,bestVX,bestVY),radMxs(bestZ,bestVX,bestVY), bestTint);
fprintf(fileID, '      - Best Rel Int Time Loc (fibz = %i, vX = %i, vY = %i)\n', bestZ, bestVX, bestVY);
fclose(fileID);

fprintf('\nFull rel. tint. matrix\n')
relTints = permute(relTints, [3 2 1]);
disp(relTints);

% (OPTIONAL) print the min and max in each frame
% fprintf('\nFull mins matrix\n')
% disp(permute(radMns,[3,2,1]));
% 
% fprintf('\nFull maxs matrix\n')
% disp(permute(radMxs,[3,2,1]))

end
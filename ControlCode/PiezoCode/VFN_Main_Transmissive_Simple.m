%{ 
 Main_Transmissive_Simple
 Code for doing a normal scan with the VFN Transmissive bench
 - Uses the Piezos (MDT) for fiber X and Y 
 - Uses Zabers for fiber Z, vortex X, and vortex Y
 - Uses the redPM and a Zaber for in-beam normalization
 - Uses the Femto for through-fiber power measurements

 * Allows normalization using the redPM or an MMF and the Femto
   - Also possible to normalize at all
   *** The MMF values reported in the fits files are now -9999 if MMF is
        not used. This follows the same idea as the redPM when not used.
   * Code also does not share pmRead1 between MMF and pmRead anymore
 * Outputs cubes of the same dimensionality as the PiezoScanZabFoc code
   - ie: [X, Y, Nread, Focus, VortX, VortY]
 * Uses functions from VFN Analysis Library to interpret the data immediately
   - Like this, the values reported by this code are the same as well see when
    running the analysis code later (assuming same analysis params are used)
   - Full analysis (Fresnel losses, prop. losses, lens losses, etc.) is done.
 * Automatically prints values and results into the DataNotes file 
 * Uses the new control libraries for the VFN devices
 * Calculates radially-averaged profile to determine proper planet coupling
 * Calculates rel. int. time using the null and the rad. avg. planet coupling 
   - Will optionally iterate through all frames and find the best rel. int.
 * Tries to not move the Zabers unless absolutely necessary.

 Author: Daniel Echeverri (dechever@caltech.edu)
 Origin Date: 03/02/2020

%}

%% add the necessary functions
addpath(genpath('C:\Users\Daniel Echeverri\Documents\MATLAB\VFN-Lab\ControlCode'));
addpath(genpath('C:\Users\Daniel Echeverri\Documents\MATLAB\VFN-Lab\AnalysisCode\AnalysisLib'));
addpath(genpath('C:\Users\Daniel Echeverri\Documents\MATLAB\VFN-Simulations\VFNlib'));

close all
clear all

%% General settings 

% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\200701_COV3';

%-- Experiment name for figures
expNm = '10Don_785LasDeepScan2';

%-- Custom auto save text
% Change at every run
run_msg = 'Check first decent spot in #9.';
% Change when doing new set of experiments
set_msg = 'Check bench performance after realignment.';
las_msg = 'Laser 785nm, 7.0mW.';   % Laser/Varia settings
flt_msg = 'N/A';                     % Optical Filter used
vtx_msg = 'JPL Poly 550nm Charge 2';            % Vortex in use

%-- Final analysis stuff
isRadAvgAnalysis = true;    % calc rel. int. time in all frames using rad avg

%~~ ZABER STUFF 
% Flag to enable vortex scan.
isVortScan = false;

% Distance to move zabers for backlash removal [in mm]
zabBacklash= 0.015;         %Note: affects vortex and fiber z-axis

% Vortex center of scan [in mm]
VXcenter = 03.312000;       % 
VYcenter = 13.917000;       %!! must be >9.9

% Vortex scan properties
%To include center values as a point in scan; use ODD number of points
VXpoints = 3;           % Number of X points in scan
VYpoints = VXpoints;    % Number of Y points in scan

% Vortex step params in [mm]  
vStepRange = min(0.008,0.8);   %Vortex will be scanned +/- this value
                    % ie. scan will be: 
                    % [VXcen-vStepRa to VXcen+vStepRa] w/ Vxpoi steps
                    
% Fiber Focus Scan properties
%To include center value as a point; use ODD number of points
Zcenter = 7.032000;     % Focus (Z) center point. !! must be < 9
Zpoints = 1;            % Number of focci taken (exact number; no longer +1)
ZStepSize = 0.007;      % Step size in mm
%~~ END ZABER STUFF 

%~~ PIEZO STUFF 
% Fiber center of scan in Volts
Xcenter = 75.00;
Ycenter = 100.00;

% Fiber scan properties
%To include center values as a point in scan; use EVEN number of points
Xpoints = 50;% number of X points in scan (actual scan will have +1)
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle

% Fiber step sizes in Volts
refStep   = 10; % refStep is the step size for a full PSF with a 10*10 grid
StepSize  = refStep/(Xpoints/10);
backlash  = 10; % Kept backlash in to see if it did anything
%~~ END PIEZO STUFF 

%~~ FEMTO POWER METER STUFF 
isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
Nread   = 100;      % power samples at a given locaiton
Nrate   = 1000;     % rate of scan [samples/second]
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
nrmValReport = 1/44.1;   % @635nm typ = 14.08; @780nm typ = 44.1
                        % Valu for fibertip power when isPMNorm=false. Only
                        % affects the printed values at end; NORM'd cube
                        % values are still un-normed. SET =1 to not norm
                        % the printed values.
pmCalWvl = 785; % Wavelength for redPM for normalization
%~~ END RED (NORM) POWER METER STUFF

%% Zaber Setup
fprintf('---Performing Zaber Setup\n')
if isVortScan || (Zpoints > 1) || isPMNorm
    % Enable zabers since they will be used at some point
    isZab = true;
end

%~~ Connect to Zabers to query position
VFN_setUpBenchZabers; % Instantiate and rename the zabers

%- Make sure pmX is out of the way
if isfield(Zabs, 'pmX')
    VFN_Zab_move(Zabs.pmX, 25);
elseif isPMNorm
    error('PM Norm is requested but no pmX Zaber exists.')
end

%- Query positions
if (VFN_Zab_getPos(Zabs.vortX) ~= VXcenter) || (VFN_Zab_getPos(Zabs.vortY) ~= VYcenter)
    % vort not at cent. pos; if not doing vort scan, enable zabs to correct pos later
    isZab = true;
end
if VFN_Zab_getPos(Zabs.fibZ) ~= Zcenter
    % foc not at cent. pos; if not doing foc scan, enable zabs to correct pos later
    isZab = true;
end

%- Disable Zabers for rest of code if not needed anymore
if ~isZab 
    VFN_cleanUpZabers;
end
    

%% MDT (piezo) Setup
fprintf('---Performing MDT Setup\n')
% Connect to the MDT Piezos
VFN_setUpMDT;

% Query MDT on each pertinent axis to ensure correct reads later
    % It seems the very first read after connection might be incorrect
VFN_MDT_getPos(MDT, 'x')
VFN_MDT_getPos(MDT, 'y')

%% Femto setup
fprintf('---Performing Femto Setup\n')
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

% Set FMTO gain to user-provided value (known start value for code)
VFN_FMTO_LUCI_setGain(FMTO_scale);

fprintf('Current Gain setting: %i\n', FMTO_scale)

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
datName = '\DataNotes.txt';
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
fprintf(fileID, '        scan center    (%5.2f V,%6.2f V, %8.6f mm)\r\n', trash);
trash = [Xpoints; Zpoints; refStep; refStep; Xpoints; ZStepSize];
fprintf(fileID, '        Scan params:   xpoi = %d, zpoi = %d, refStep = %d, stepsize = %d/(%d/10); zstep = %4.2f mm\r\n', trash);
fprintf(fileID, '        isAutoScale:   %s\r\n', mat2str(isAutoScale));
trash = [FMTO_scale; Nread; Nrate];
fprintf(fileID, '        Femto:         FMTO_scale = %d; Nread = %d; Nrate = %d\r\n', trash);
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
meas= nan(length(distX),length(distY),Nread,length(distZ),length(distVX), length(distVY)); 
measScl = meas;         %Matrix for semi-processed data
sclSz   = size(meas); sclSz(3) = 3;
scales  = nan(sclSz);   %Matrix for scales, gains, and biases

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
      % check for collisions first
    if (distVY(1)-zabBacklash) >= 9.9
        % VortY axis will not collide in the next move
        if isZab && (distVY(1) ~= VFN_Zab_getPos(Zabs.vortY))
            % Remove y backlash now that motion is confirmed safe
            VFN_Zab_move(Zabs.vortY, distVY(1)-zabBacklash);
        end
    else
        error('vortY position will collide with pmX')
    end
    
    for b = 1:length(distVY)
        if isZab && (distVY(b) ~= VFN_Zab_getPos(Zabs.vortY))
            %move vortex to new Y position
            fprintf('  VortY pos: %f\n',VFN_Zab_move(Zabs.vortY, distVY(b)));
        end

        % Remove zaber focus (Z) backlash
          % check for collisions first
        if (distZ(1)-zabBacklash) <= 9
            % fibZ axis will not collide in the next move
            if isZab && (distZ(1) ~= VFN_Zab_getPos(Zabs.fibZ))
                % Remove z backlash now that motion is confirmed safe
                VFN_Zab_move(Zabs.fibZ, distZ(1)-zabBacklash);
            end
        else
            error('fibZ position will collide with lens mount')
        end
        
        for k=1:length(distZ)
            if isZab && (distZ(k) ~= VFN_Zab_getPos(Zabs.fibZ))
                % Move fiber in focus
                VFN_Zab_move(Zabs.fibZ, distZ(k));
            end
            
            % Remove MDT (X) backlash
            VFN_MDT_move(MDT, 'x', distX(1)-backlash);
            
            for j=1:Xpoints+1
                % Print current scan status
                fprintf('Vort %i of %i | Foc %i of %i | Column %i of %i\n', ...
                        (a-1)*length(distVY)+b,length(distVY)*length(distVX),...
                        k, length(distZ), j, length(distX))
                
                % Move MDT(X)
                VFN_MDT_move(MDT, 'x', distX(j));
                
                % Remove MDT (Y) backlash
                VFN_MDT_move(MDT, 'y', distY(1)-backlash);
                
                for i=1:Ypoints+1
                    % Move MDT (Y) 
                    VFN_MDT_move(MDT, 'y', distY(i));
                    
                    % Let femto voltage settle at new position
                    pause(Delay);
                                       
                    % Take a sample reading at current power
                    read    = startForeground(s);
                    
                    % Perform auto scaling if desired
                    if isAutoScale
                        % Save old FMTO_scale for comparison after autoGain
                        old_scale = FMTO_scale;
                        % Modify gain accordingly
                        [FMTO_scale, s] = VFN_FMTO_LUCI_setAutoGain(s, mean(read), FMTO_scale);
                        % Check if re-read is needed
                        if old_scale ~= FMTO_scale
                            % FMTO_scale changed so the gain changed
                            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
                            read    = startForeground(s);
                            if FMTO_scale > 9
                                warning('Gain >9: %i',FMTO_scale)
                            end
                        end
                    end
                    if find(read>10)
                        warning('Power is too high')
                    end
                    
                    % Store data
                    meas(j,i,:,k,a,b) = read;
                    scales(j,i,1,k,a,b) = FMTO_scale;
                    scales(j,i,2,k,a,b) = gnFact^-(FMTO_scale-6);
                    
                    % Get Bias at current gain setting
                    locBias = VFN_FMTO_getBias(FMTO_scale);
                    
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
    mmf_read    = startForeground(s);
    % Perform auto scaling if desired
    if isAutoScale
        % Save old FMTO_scale for comparison after autoGain
        old_scale = FMTO_scale;
        % Modify gain accordingly
        [FMTO_scale, s] = VFN_FMTO_LUCI_setAutoGain(s, mean(mmf_read), FMTO_scale);
        % Check if re-read is needed
        if old_scale ~= FMTO_scale
            % FMTO_scale changed so the gain changed
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
            mmf_read    = startForeground(s);
            if FMTO_scale > 9
                warning('Gain >9: %i',FMTO_scale)
            end
        end
    end
    if find(mmf_read>10)
        warning('Power is too high')
    end
    
    % Record femto gain setting during MMF read
    MMF_FMTO_scale = FMTO_scale;    
else
    %-- Set mmf_read and MMF_FMTO_scale to -9999 if MMF is not used
    mmf_read = ones(Nread,1)*-9999;
    MMF_FMTO_scale = -9999;
end

%-- Use the red power meter to measure the power
if isPMNorm
    %-- Move redPM into the beam
    % Check that no collisions will occur from motion
    if VFN_Zab_getPos(Zabs.vortY) < 9.9
        error('Vortex Y-Axis is in the way')
    end
    %Move if no collision
    VFN_Zab_move(Zabs.pmX, 0);

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
else
    %-- Set pmRead to -9999 if PM is not to be read
    pmRead = ones(pmNread,1)*-9999;
end

%% close the connections to all devices
% Set Femto back to 10^6 
FMTO_scale = 6;
VFN_FMTO_LUCI_setGain(FMTO_scale);
fprintf('\n Femto Gain set to %i',FMTO_scale);
VFN_cleanUpFMTO;

%Set MDT axes to central position
VFN_MDT_move(MDT, 'x', Xcenter-backlash);
VFN_MDT_move(MDT, 'x', Xcenter);
VFN_MDT_move(MDT, 'y', Ycenter-backlash); 
VFN_MDT_move(MDT, 'y', Ycenter); 
% Close connection to the MDT
VFN_cleanUpMDT

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
    
    %-- Move fiber back to Zcenter
    if Zcenter ~= VFN_Zab_getPos(Zabs.fibZ)
        VFN_Zab_move(Zabs.fibZ, Zcenter-zabBacklash);
        % Center Zab (Z)
        VFN_Zab_move(Zabs.fibZ, Zcenter);
    end
    
    % Print position of vortex
    fprintf('\nVortX pos: %f  VortY pos: %f  \n    fibZ pos: %f',...
        VFN_Zab_getPos(Zabs.vortX), ...
        VFN_Zab_getPos(Zabs.vortY), ...
        VFN_Zab_getPos(Zabs.fibZ));
    
    
    %-- Move Calibration PM out of the beam
    if isPMNorm
        fprintf('  |  PM Zab pos: %f\n', VFN_Zab_move(Zabs.pmX, 25));
    else
        fprintf('\n');  % terminate line from zaber recentering
    end
    
    %-- Clean up
    VFN_cleanUpZabers;
end

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

%-- Define a separate matrix which normalizes the data by power on the 
if isPMNorm
    %-- When actual PM reading was done, normalize as usual
    % Account for loss from fiber lens
        % 3.4% reflection spec at 780nm, 0.3% reflection measured at 635nm
    pmRead1 = mean(pmRead)*(1-0.034);%0.997;    % Updated on 2/8/19 based on 10/26/18
    % Convert the calibration PM value from W to uW
    pmRead1 = pmRead1*10^6;
    % Translate the calibration PM value from uW to V (at 10^6)
    fmtoSnsScl = VFN_getFmtoResponsivity(pmCalWvl);
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
    scl2 = gnFact^-(MMF_FMTO_scale-6);
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
            imagesc((distX-Xcenter),(distY-Ycenter),transpose(measScl2(:,:,:,k,a,b)));
            titStr = sprintf('PWR Thru Fib');
            if isZab
                % Append zaber positions if moved
                titStr = [titStr sprintf(', Z= %0.5f mm, vort=(%0.5f, %0.5f)',distZ(k), distVX(a),distVY(b))]; %#ok<AGROW,*UNRCH>
            end
            title(titStr)
            axis xy;axis image%axis square;
            cbr = colorbar;
            ylabel(cbr, 'Normed Coupling')
            set(gca,'TickDir','out')
            xlabel('Fiber X position [V]')
            ylabel('Fiber Y position [V]')
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

% % Optional plot of line profile across center row
% figure()
% plot(measScl2(:,round(Ypoints/2+1)))
% xlabel('horizontal position of the fiber in the  focal plane in micrommeter')
% ylabel('power in mWatts')
% title(strcat('scan horizontal at ',num2str(Delay),'sec of delay'))
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
NormCubeNm    = [svFld filesep expNm '_' tag '_NORMREDU.fits'];

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
fits.writeKey(fitmap, 'Nread ',  Nread);
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
fits.writeKey(fitmap, 'isAutScl', logical(isAutoScale));
ind = strfind(GnCubeNm ,'\');
fits.writeKey(fitmap, 'DatCube ', GnCubeNm(ind(end-1)+1:end));
ind = strfind(SRCubeNm ,'\');
fits.writeKey(fitmap, 'SRCube ', SRCubeNm(ind(end-1)+1:end));
ind = strfind(NormCubeNm ,'\');
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
fits.writeKey(fitmap, 'CNTRLCD', 'Main_TransSimple');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(meas,permOrd));
fits.closeFile(fitmap);

%-- SAVE GAIN SETTINGS
fitmap= fits.createFile(GnCubeNm);
fits.createImg(fitmap,'double',size(permute(scales, permOrd)));
%header data
ind = strfind(datCubeNm,'\');
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
ind = strfind(datCubeNm,'\');
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
ind = strfind(datCubeNm,'\');
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
fprintf('\nMin in center region:    %f',mnVal)
fprintf(['\n Min Location: ', ...
         '\n  X index = %3i, X   = %6.2f     V', ...
         '\n  Y index = %3i, Y   = %6.2f     V', ...
         '\n  Foc ind = %3i, Foc =   %0.5f  mm', ...
         '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
         '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
         mnInd(1), distX(mnInd(1)), mnInd(2), distY(mnInd(2)), ...
         mnInd(4), distZ(mnInd(4)), mnInd(5), distVX(mnInd(5)), ...
         mnInd(6), distVY(mnInd(6)));


%-- Find overall max
[mxVal2, mxInd2] = VFN_An_getEta_p(measScl);
fprintf('\nOverall Max:             %f',mxVal2)
fprintf('\n Overall Max Location: (%i, %i, %i, %i, %i, %i)', ...
        mxInd2(1), mxInd2(2), mxInd2(3), mxInd2(4), mxInd2(5), mxInd2(6));
    
%-- Find max in frame with min value
meas3 = measScl(:,:, mnInd(3), mnInd(4), mnInd(5), mnInd(6)); %Frame with min
mxVal = max(meas3(:)); 
fprintf('\nMax in min Frame:        %f', mxVal)
fprintf('\nMax/Min, min Frame: %f', mxVal/mnVal)
if ~(isPMNorm || isMMFNorm)
    fprintf('\nWARN: redPM not used; fib tip power provided by user');
end

%-- Get normed values (used measScl to get only gain-corrected values above)
eta_s = mnVal*nrmVal;
eta_p = mxVal*nrmVal;

fprintf('\nPower on Fiber tip:      %f', 1/nrmVal);
fprintf('\nMin, normed (eta_s):  %f', eta_s);
fprintf('\nMax, normed (eta_p):  %f\n', eta_p);

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
        '\n  X index = %3i, X   = %6.2f     V', ...
        '\n  Y index = %3i, Y   = %6.2f     V', ...
        '\n  Foc ind = %3i, Foc =   %0.5f  mm', ...
        '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
        '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
        mxInd2(1), distX(mxInd2(1)), mxInd2(2), distY(mxInd2(2)), ...
        mxInd2(4), distZ(mxInd2(4)), mxInd2(5), distVX(mxInd2(5)), mxInd2(6), distVY(mxInd2(6)));
    
    fprintf('\n\n- Max:                        %f', mxVal2)
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

end
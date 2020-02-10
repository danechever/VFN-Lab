%{ 
 Piezo Scan
 Code for scanning the fiber in the focal plane using the Piezo actuators
 It implements some of the things learned in Nicolas's experiments

 Removed all movement limit stuff since the Pzt's can only move 14um so
    they are extemely unlikely to collide with anything

*** This version saves a [X, Y, Nread, Focus] dimensional cube in the fits
      file instead of saving just the mean of the Nreads. 

*** This version lets you control the vortex position using the zabers

*** This version enables vortex raster scans

*** This version automatically changes the gain on the femto to give a high
      dynamic range on the reads.

*** This version uses the red thorlabs power meter to normalize the values

*** [04/17/19] Modified to make red PM read for normalization optional
    - Added isPMNorm flag to choose whether red PM is read or not
    - When red PM is not read, -9999 is saved as pmRead value in the data
        cube and the normed cube is "normed" by 1. Also, all values
        reported at end of script are normed by user provided value.
    - Added isPMNorm flag to keywords of data and normed cubes
    * Also, fixed path in first line to reflect the fact that all the
        control code is now in a git repo

*** [05/15/19] This version uses a zaber for the fiber focus (Z) axis

%}
% Programm to run the piezo's in the VFN experiment to see the donut

%% add the necessary functions
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode'));

close all
clear all

%% General settings 

% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\020620_FLE7';

% Experiment name for figures. If the first number is single digit,
expNm = '6Don_780VarVortScan2'; % please pad with zero, ex. 01Don_
message = 'Check full donut coupling before another vortex scan.';

%~~ ZABER STUFF 
% Flag to enable zaber motion (including centering)
isZab = true;   % Should be true if fiber focus scan is desired
isVortScan = false;

% Vortex center of scan properties [in mm] %11.963700 in, 17.485000 out
VXcenter = 12.016700;      % Value when PM stage was added: 12.271445mm
VYcenter = 12.895000;  %!! must be >12.334083
%VYcenter 12.880000 in, 17.500000 out
Vbacklash= 0.005;          % Note: this is also used for fiber z-axis

% Vortex scan properties
%To include center values as a point in scan; use ODD number of points
% !!! THIS IS VY !!!
VYpoints = 3;% Number of Y points in scan
% !!! THIS IS VX !!!
VXpoints = VYpoints; % Number of X points in scan

% Vortex step params in [mm]  
vStepRange = min(0.015,0.8);   %Vortex will be scanned +/- this value
                    % ie. scan will be: 
                    % [VXcen-vStepRa to VXcen+vStepRa] w/ Vxpoi steps

% Fiber Focus Scan properties
Zcenter = 1.743000;     % Focus (Z) center point %% Change by BC from 3.008700
Zpoints = 1;            % Number of focci taken
ZStepSize = 0.007;      % Step size in mm
%~~ END ZABER STUFF 

%~~ PIEZO STUFF 
% Fiber center of scan in Volts
Xcenter = 80.00; % clockwise moves the center to the left
Ycenter = 81.00; % clockwise moves the center down

% Fiber scan properties
%To include center values as a point in scan; use EVEN number of points
Xpoints = 50;% number of X points in scan (actual scan will have +1)
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle

% Fiber step sizes in Volts
refStep   = 12; % refStep is the step size for a full PSF with a 10*10 grid
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
delLookUp   = [[5 1]; [1 0.5]; [0 0.3]];  % Look-up table for delay values
delInd      = find(delLookUp<=StepSize,1);  % index of delay to use
Delay       =delLookUp(delInd,2); 
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)
%~~ END FEMTO POWER METER STUFF 

%~~ FEMTO NORMALIZATION STUFF
isMMFNorm = false;
%~~ END FEMTO NORMALIZATION

%~~ RED (NORM) POWER METER STUFF 
pmNread = 100;      % Number of samples to take
isPMNorm = true;   % Flag to mark whether this PM should be read
                        % This flag is useful in case the PM is not in the
                        % system. -9999 will replace the pmRead values.
pmNormReport = 44.1;   % @635nm typ = 14.08; @780nm typ = 44.1
                        % Valu for fibertip power when isPMNorm=false. Only
                        % affects the printed values at end; NORM'd cube
                        % values are still un-normed. SET =1 to not norm
                       % the printed values.
pmCalWvl = 780; % Wavelength for redPM for normalization
                        
%~~ END RED (NORM) POWER METER STUFF
%% DataNotes automation pt1
% Write the filename
datName = strcat('\DataNotes.txt');
datFl = strcat(svFld, datName);
fileID = fopen(datFl, 'a');
fprintf(fileID, '- %s\r\n', expNm);

% General comment
 fprintf(fileID, '    %s\r\n', message);
fprintf(fileID, '    * Trying to improve planet coupling on varia.\r\n');

% Scan Params section
fprintf(fileID, '    Scan params:\r\n');
fprintf(fileID, '        isZab:         %s\r\n', mat2str(isZab));
fprintf(fileID, '        isVortScan:    %s\r\n', mat2str(isVortScan));%    
trash = [VXcenter; VYcenter];
fprintf(fileID, '        Vortex pos:    x = %9.6fmm    y = %9.6fmm\r\n', trash);
fprintf(fileID, '        Vbacklash:     %5.3fmm\r\n', Vbacklash); 
trash = [VYpoints; VXpoints; vStepRange];
fprintf(fileID, '        v scan param:  vYpoi = %d, vXpoi = %d, vStepRange = %4.2f\r\n', trash);
fprintf(fileID, '        Laser power:   Varia 780nm, 80%% power, 3nm BW.\r\n');
fprintf(fileID, '        Optical Filt:  600nm Longpass\r\n');
fprintf(fileID, '        Vortex Mask:   JPL Poly 550nm Charge 2\r\n');
trash = [Xcenter; Ycenter; Zcenter];
fprintf(fileID, '        scan center    (%5.2f V,%6.2f V, %8.6f mm)\r\n', trash);
trash = [Xpoints; Zpoints; refStep; refStep; Xpoints; ZStepSize];
fprintf(fileID, '        Scan params:   xpoi = %d, zpoi = %d, refStep = %d, stepsize = %d/(%d/10); zstep = %4.2f mm\r\n', trash);
fprintf(fileID, '        isAutoScale:   %s\r\n', mat2str(isAutoScale));
trash = [FMTO_scale; Nread; Nrate];
fprintf(fileID, '        Femto:         FMTO_scale = %d; Nread = %d; Nrate = %d\r\n', trash);
trash = [pmNormReport; pmCalWvl];
fprintf(fileID, '        redPM norm:    isPMNorm = %s; pmNormReport = %5.2f; pmCalWavelength = %dnm\r\n', mat2str(isPMNorm), trash);
fprintf(fileID, '        isMMFNorm: %s\r\n\r\n', mat2str(isMMFNorm));
fclose(fileID);

%~~ END DataNotes pt1 automation
%% Zaber Setup
if isVortScan && ~isZab
    error('isVortScan is true but vortex motion is disabled with isZab')
end

if (Zpoints > 1) && ~isZab
    error('Fiber scan is desired but zaber motion is disabled with isZab')
end

if isZab
    VFN_setUpZabers; % Instantiate the zabers

    % Relabel the variables for clarity
    if exist('ax14', 'var')
        vortX = ax14;
        clear ax14
    end
    if exist('ax63', 'var')
        vortY = ax63;
        clear ax63
    end
    if exist('ax93', 'var')
        fibZ = ax93;
        clear ax93
    end
    if exist('ax214', 'var')
        pmX = ax214;
        clear ax214
    end
    
    if isPMNorm
        % Move red PM out of beam
        VFN_Zab_move(pmX, 25);
    end
    
end

%% PI Setup
%Create serial port and define communication properties
    %COM# changes between computers; must check in Device Manager
VFN_setUpPIStage

%****NEED TO CHECK VALUES FROM QFRF ARE 0/1******
%Check to make sure axes are referenced.
for idx = 1:length(axis)
    if PIDevice.qFRF(char(axis(idx))) == 0
        temp=input("Axis "+idx+" needs to reference. Ok to move? y/n [n]", "s");
        if temp=="n"
            VFN_cleanUpPIStage;
            error("User unwilling to reference stage")
        elseif temp=="y"
            VFN_PIStage_home(idx, PIDevice)
        end
    end
end

%% Femto setup
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

%if isAutoScale
%    % Take a sample reading at current power
%    readVal = mean(startForeground(s));
%    % Modify gain accordingly
%    [FMTO_scale, s] = VFN_FMTO_setAutoGain(s, readVal, FMTO_scale);
%else
    VFN_FMTO_LUCI_setGain(FMTO_scale);
%end

fprintf('Current Gain setting: %i\n', FMTO_scale)

%% Red Thorlabs PM setup

%-- Connect to red PM only if it will be used
if isPMNorm
    % Find a VISA-USB object.
        % NOTE: check the SN on the device. The 2 student-lab PM100D's are:
        %   P0015519    and     P0015560
        % our VFN PM100D is: P0022602
    obj1 = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x1313::0x8078::P0022602::0::INSTR', 'Tag', '');

    % Create the VISA-USB object if it does not exist
    % otherwise use the object that was found.
    if isempty(obj1)
        % NOTE: check the SN on the device. The 2 student-lab PM100D's are:
        %   P0015519    and     P0015560
        % our VFN PM100D is: P0022602
        obj1 = visa('NI', 'USB0::0x1313::0x8078::P0022602::0::INSTR');
    else
        fclose(obj1);
        obj1 = obj1(1);
    end

    % Connect to instrument object, obj1.
    fopen(obj1);
    
    % Set operating wavelength for calibration
    fprintf(obj1, ['sense:correction:wavelength ' num2str(pmCalWvl)]);
end
%% Perform scan
distX=(Xcenter-StepSize*Xpoints/2):StepSize:(Xcenter+StepSize*Xpoints/2);
distY=(Ycenter-StepSize*Ypoints/2):StepSize:(Ycenter+StepSize*Ypoints/2);
if Zpoints ~= 1
    distZ=(Zcenter-ZStepSize*Zpoints/2):ZStepSize:(Zcenter+ZStepSize*Zpoints/2); 
else
    distZ=Zcenter;
end

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

%Preallocate data matrices
meas=zeros(length(distX),length(distY),Nread,length(distZ),length(distVX), length(distVY)); 
measScl = meas;         %Matrix for semi-processed data
sclSz   = size(meas); sclSz(3) = 3;
scales  = nan(sclSz);   %Matrix for scales, gains, and biases
%ratio   = 1;
%fstMsFlg= true;
% measPos = zeros(Xpoints,Ypoints,Zpoints,3);

if isZab
	%remove zaber X backlash
    VFN_Zab_move(vortX, distVX(1)-Vbacklash);
end
for a = 1:length(distVX)
    if isZab
        %move vortex to new X position
        fprintf('\nVortX pos: %f',VFN_Zab_move(vortX, distVX(a)));
        %remove zaber y backlash
        if (distVY(1)-Vbacklash) >= 12.33
            VFN_Zab_move(vortY, distVY(1)-Vbacklash);
        else
            error('vortY position will collide with lens mount')
        end
    end
    for b = 1:length(distVY)
        if isZab
            %move vortex to new Y position
            fprintf('  VortY pos: %f\n',VFN_Zab_move(vortY, distVY(b)));
            % Remove zaber focus (Z) backlash
            if (distZ(1)-Vbacklash) <= 7.8
                VFN_PIStage_move('z', PIdevice, distZ(1));
            else
                error('fibZ position will collide with lens mount')
            end
        end
        for k=1:length(distZ)
            
            if isZab
                % Move fiber in focus
                VFN_PIStage_move('z', PIdevice, distZ(k))
            end
            
            VFN_PIStage_move('x', PIdevice, distX(1));
            for j=1:Xpoints+1
                fprintf('Vort %i of %i | Foc %i of %i | Column %i of %i\n', ...
                        (a-1)*length(distVY)+b,length(distVY)*length(distVX),...
                        k, length(distZ), j, length(distX))
                
                VFN_PIStage_move('x', PIdevice, distX(j))   
                
                VFN_PIStage_move('y', PIdevice, distY(1))   
                for i=1:Ypoints+1
                    
                    VFN_PIStage_move('y', PIdevice, distY(i))   
                    %read=zeros(Nread,1);
                    %readPos = zeros(Nread, 3); % nReadx3 matrix for holding zaber positions
                    %col: (x, y, z), rows (reads)
                    %tic;
                    pause(Delay);
                    %if fstMsFlg
                    %    if isAutoScale
                    %        % Take a sample reading at current power
                    %        readVal = mean(startForeground(s));
                    %        % Modify gain accordingly
                    %        [FMTO_scale, s] = VFN_FMTO_setAutoGain(s, readVal, FMTO_scale);
                    %    end
                    %    fstMsFlg = false;
                    %end
                    
                    % Take a sample reading at current power
                    read    = startForeground(s);
                    if isAutoScale
                        % Save old FMTO_scale for comparison after autoGain
                        old_scale = FMTO_scale;
                        % Save old read value for comparison after autoGain
                        old_read = mean(read);
                        % Modify gain accordingly
                        [FMTO_scale, s] = VFN_FMTO_LUCI_setAutoGain(s, old_read, FMTO_scale);
                        % Check if re-read is needed
                        if old_scale ~= FMTO_scale
                            % FMTO_scale changed so the gain changed
                            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
                            read    = startForeground(s);
                            %ratio   = ratio * (old_read/mean(read));
                            if FMTO_scale > 9
                                warning('Gain >9: %i',FMTO_scale)
                            end
                        end
                    end
                    if find(read>10)
                        warning('Power is too high')
                    end
                    %for n=1:Nread
                    %    read(n)=s.inputSingleScan; % taking power
                        %if read(n) >10
                        %    warning('Power exceeds Femto limit of 10')
                        %end
                        %readPos(n, 1) = ax14.getposition;
                        %readPos(n, 2) = ax63.getposition;
                        %readPos(n, 3) = ax93.getposition;
                        
                    %end
                    %             if ~ isempty(find(readPos(2:end,1) ~= readPos(1,1),1))
                    %                 fprintf('X pos not constant\n')
                    %             end
                    %             if ~ isempty(find(readPos(2:end,2) ~= readPos(1,2),1))
                    %                 fprintf('Y pos not constant\n')
                    %             end
                    %             if ~ isempty(find(readPos(2:end,3) ~= readPos(1,3),1))
                    %                 fprintf('Z pos not constant\n')
                    %             end
                    %             mreadPos = mean(readPos,1);
                    %test=toc;
                    meas(j,i,:,k,a,b) = read;
                    scales(j,i,1,k,a,b) = FMTO_scale;
                    scales(j,i,2,k,a,b) = gnFact^-(FMTO_scale-6);
                    switch FMTO_scale
                        case 5
                            locBias = 0.004905;%0.0033887;
                        case 6
                            locBias = 0.005382;%0.0042960;
                        case 7
                            locBias = 0.005405;%0.0044690;
                        case 8
                            locBias = 0.005350;%0.0047430;
                        case 9
                            locBias = 0.004941;%0.0068927;
                        case 10
                            locBias = -0.001120;
                        case 11
                            locBias = -0.059031;
                        otherwise
                            warning('No bias forFMTO_scale = %i\n',FMTO_scale)
                            locBias = nan;
                    end
                    scales(j,i,3,k,a,b) = locBias;
                    measScl(j,i,:,k,a,b) = read - locBias; %Subtract bias for semi-processed matrix
                    %meas(j,i,k)=mean(read);
                    %measPos(j,i,k,:) = mreadPos(:);
                end
                
            end
        end
    end
end
%meas=meas/pm_scale;

%% Measure power at Fiber Tip
%-- Use the red power meter to measure the power
% Move the PM into the beam
if isZab && isPMNorm
    % Check that no collisions will occur from motion
    if VFN_Zab_getPos(vortY) < 12.35
        error('Vortex Y-Axis is in the way')
    end
    %Move if no collision
    VFN_Zab_move(pmX, 0);
end

if isPMNorm
    %-- Actually read from red PM when flag is set. 
    % Wait for PM value to settle
    pause(2);

    % Take measurements of the power
    pmRead = nan(pmNread,1);
    for ii = 1:pmNread
        pmRead(ii)=str2num(query(obj1, 'measure:power?'));
    end
    % Account for ~0.02 loss from fiber lens
        % 3.4% reflection spec at 780nm, 0.3% reflection measured at 635nm
    pmRead1 = mean(pmRead(:))*(1-0.034);%0.997;    % Updated on 2/8/19 based on 10/26/18
else
    %-- Set pmRead to -9999 if PM is not to be read
    pmRead = ones(pmNread,1)*-9999;
end

if isMMFNorm
    disp('Switch fibers. Press any key to continue')
    %%%% When this becomes automated, this will become a PI command
    %%%% rather than a manual switch CHANGE_MARK
    pause;
    %Now the MMF should be in the beam and can collect all the light
    %Use the femto to collect the total power, same as collecting data
    mmf_read    = startForeground(s);
    if isAutoScale
        % Save old FMTO_scale for comparison after autoGain
        old_scale = FMTO_scale;
        % Save old read value for comparison after autoGain
        old_read = mean(mmf_read);
        % Modify gain accordingly
        [FMTO_scale, s] = VFN_FMTO_LUCI_setAutoGain(s, old_read, FMTO_scale);
        % Check if re-read is needed
        if old_scale ~= FMTO_scale
            % FMTO_scale changed so the gain changed
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
            mmf_read    = startForeground(s);
            %ratio   = ratio * (old_read/mean(read));
            if FMTO_scale > 9
                warning('Gain >9: %i',FMTO_scale)
            end
        end
    end
    if find(mmf_read>10)
        warning('Power is too high')
    end
    
    %store the relevant values in the same manner as before
    scl1 = FMTO_scale;
    scl2 = gnFact^-(FMTO_scale-6);
    switch FMTO_scale
        case 5
            locBias = 0.004905;%0.0033887;
        case 6
            locBias = 0.005382;%0.0042960;
        case 7
            locBias = 0.005405;%0.0044690;
        case 8
            locBias = 0.005350;%0.0047430;
        case 9
            locBias = 0.004941;%0.0068927;
        case 10
            locBias = -0.001120;
        case 11
            locBias = -0.059031;
        otherwise
            warning('No bias forFMTO_scale = %i\n',FMTO_scale)
            locBias = nan;
    end
    scl3 = locBias;
    mmf_norm = mmf_read-locBias;
    pmRead = mmf_norm;
end


%% close the connections to all devices
% Set Femto back to 10^6 
FMTO_scale = 6;
VFN_FMTO_LUCI_setGain(FMTO_scale);
fprintf('\n Femto Gain set to %i',FMTO_scale);
VFN_cleanUpFMTO;

%Set all three axes to starting position
VFN_PIStage_move('x', PIdevice, Xcenter)
VFN_PIStage_move('y', PIdevice, Ycenter)
VFN_PIStage_move('z', PIdevice, Zcenter)

%clean up PI stage
VFN_cleanUpPIStage

if isZab
    
    if isPMNorm
        %-- Move Calibration PM out of the beam
        fprintf('  |  PM Zab pos: %f\n', VFN_Zab_move(pmX, 25));
    else
        fprintf('\n');  % terminate line from vortex recentering
        
    %-- Return vortex to center position
    VFN_Zab_move(vortX, VXcenter-Vbacklash);
    fprintf('\nVortX pos: %f',VFN_Zab_move(vortX, VXcenter));
    VFN_Zab_move(vortY, VYcenter-Vbacklash);
    fprintf('  Vort Y pos: %f',VFN_Zab_move(vortY, VYcenter));
    %-- Move fiber back to Zcenter
    VFN_PIStage_move('z', PIdevice, Zcenter)
    end
        
    %-- Clean up
    VFN_cleanUpZabers;
    if exist('vortX', 'var')
        clear vortX
    end
    if exist('vortY', 'var')
        clear vortY
    end
    if exist('fibZ', 'var')
        clear fibZ
    end
    if exist('pmX', 'var')
        clear pmX
    end
end

%-- Disconnect from calibration PM
if isPMNorm
    %-- Only disconnect if we connected earlier
    fclose(obj1);
end

%% Scale data accordingly to compensate for varied gain
%-- This is for plotting and analysis; the raw data is still stored as well

%** Bias is already subtracted in the main for loop sectio

%-- Define scaled matrix; equal to mean of meas unless gain was varied
measScl = mean(measScl,3).*scales(:,:,2,:,:,:);

%-- Define a separate matrix which normalizes the data by power on the 
if isPMNorm
    %-- When actual PM reading was done, normalize as usual
    % Convert the calibration PM value from W to uW
    pmRead1 = pmRead1*10^6;
    % Translate the calibration PM value from uW to V (at 10^6)
    fmtoSnsScl = VFN_getFmtoResponsivity(pmCalWvl);
    pmRead1 = pmRead1*fmtoSnsScl;     % V/uW conversion is ~ 0.63 @633; ~0.9 @ 780
else
    %-- When PM was not read, set normalization value to 1
    pmRead1 = 1;
end

if isMMFNorm
    pmRead1 = mean(mmf_norm).*scl2;
    %fmtoSnsScl = VFN_getFmtoResponsivity(pmCalWvl);
    %pmRead1 = pmRead1*fmtoSnsScl; 
end
% Create the matrix and normalize by total power on fib tip
measScl2 = measScl/pmRead1;


%% Display of data
for a = 1:length(distVX)
    for b = 1:length(distVY)
        for k = 1:length(distZ)
            figure();
            colormap('gray');
            %Not sure why transpose and "axis xy" are needed for correct orient.
            imagesc((distX-Xcenter),(distY-Ycenter),transpose(measScl(:,:,:,k,a,b)));
            titStr = sprintf('PWR Thru Fib');
            if isZab
                % Append zaber positions if moved
                titStr = [titStr sprintf(', Z= %0.5f mm, vort=(%0.5f, %0.5f)',distZ(k), distVX(a),distVY(b))]; %#ok<*UNRCH>
            end
            title(titStr)
            axis xy;axis image%axis square;
            cbr = colorbar;
            ylabel(cbr, 'Power [V]')
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

%% focus test
% If you wanna find the focus, do the normal test with the following
% section uncommented. This section will create a plot of the max power
% found in each image to see where you are closer from the focus
% for k = 1:Zpoints
%     maxVals(k) = max(max(meas(:,:,k)));
% 
% end
% figure;
% plot(distZ(1:end-1),maxVals,'-o');
% xlabel('Z (micrometer)')
% ylabel('Power')
% if you wanna look for the vortex foxus, do the same thing after replacing
% max(max by min(min)). For this one to work you need to have the donut a little
% bit bigger than the field.
%min(min(meas))

%% looking for time constant
% if you wanna calibrate the delay for your powermeter you have to : find center of the PSF, put delay to 0, ask 2 point in X, 0 points in Y, Nread = 1000, adjsut the stepsize.
% then uncomment the following code :

% figure()
% plot(0:test/Nread:test-test/Nread,read*1000)
% xlabel('time in seconds')
% ylabel('power in milliwats')
% title(strcat('Power read by the program at max speed (no delay) with 1000 reading, stepsize of ',num2str(StepSize*1000),' microns'))
% % rms of the values
% rmsread=zeros(Nread/5-1,1);
% for n=1:Nread/5-1
%     rmsread(n)=rms(read(n*5:n*5+5)*1000)-mean(read(n*5:n*5+5)*1000);
% end
% figure()
% plot(0:test/Nread*5:test-2*test/Nread*5,rmsread)
% xlabel('time in seconds')
% ylabel('difference to the mean in milliwats')
% title(strcat('rms error, stepsize of ',num2str(StepSize*1000),' microns'))
figure()
plot(measScl(:,round(Ypoints/2+1)))
xlabel('horizontal position of the fiber in the  focal plane in micrommeter')
ylabel('power in mWatts')
title(strcat('scan horizontal at ',num2str(Delay),'sec of delay'))
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
fits.writeKey(fitmap, 'MMFMean', mean(pmRead));
fits.writeKey(fitmap, 'MMFSTD', std(pmRead));
fits.writeKey(fitmap, 'NormWvl', pmCalWvl);
fits.writeKey(fitmap, 'NormScl', fmtoSnsScl);
fits.writeKey(fitmap, 'CNTRLCD', 'PZScan_ZabFoc');

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
fits.writeKey(fitmap, 'NAX3', 'mean(meas-bias)*Gain/pmRead');
fits.writeKey(fitmap, 'NAX4', 'Focus');
fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');
fits.writeKey(fitmap, 'NormFlag', logical(isPMNorm));
fits.writeKey(fitmap, 'MFFlag', logical(isMMFNorm));
fits.writeKey(fitmap, 'NormMean', mean(pmRead));
fits.writeKey(fitmap, 'NormSTD', std(pmRead));
fits.writeKey(fitmap, 'MMFMean', mean(pmRead));
fits.writeKey(fitmap, 'MMFSTD', std(pmRead));

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(measScl2,permOrd));
fits.closeFile(fitmap);



%% Report results

%-- Use user-provided norm value if desired
if ~(isPMNorm || isMMFNorm)
    pmRead1 = pmNormReport;
end

% Take average of all the Nread
meas1 = measScl;
% Crop the matrix to middle (cropVal-2)/cropVal region to be in donut null for min
%eta_4a, ind4a = VFN_An_GetEta_s(meas1)
cropVal = 7;
xmin = max(floor(length(distX)/cropVal),1);
ymin = max(floor(length(distY)/cropVal),1);
xmax = (floor((cropVal-1)*length(distX)/cropVal+1));
ymax = (floor((cropVal-1)*length(distY)/cropVal+1));
meas2 = meas1(xmin:xmax,ymin:ymax,:,:,:,:);

% Take the min of the cropped region along all dimensions 
%mnVal= min(min(min(meas1(Xpoints/5:4*Xpoints/5+1,Ypoints/5:4*Ypoints/5+1,:))));
[mnVal, mnInd] = min(meas2(:)); %get min and linear index
[I1, I2, I3, I4, I5, I6] = ind2sub(size(meas2), mnInd); %linear ind to subind
fprintf('\nMin in center region:    %f',mnVal)
fprintf(['\n Min Location: ', ...
         '\n  X index = %3i, X   = %6.2f     V', ...
         '\n  Y index = %3i, Y   = %6.2f     V', ...
         '\n  Foc ind = %3i, Foc =   %0.5f  mm', ...
         '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
         '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
         I1+xmin-1, distX(I1+xmin-1), I2+ymin-1, distY(I2+ymin-1), ...
         I4, distZ(I4), I5, distVX(I5), I6, distVY(I6)); %ind41(2), (1), (4)...


%Take overall max
%mxVal = max(max(max(meas1)));
%ind = floor(find(meas1>= mxVal)/((Xpoints+1)*(Ypoints+1)));
%fprintf('\nFrame for max:        %d',ind)
[mxVal2, mxInd2] = max(meas1(:));
[mxI1, mxI2, mxI3, mxI4, mxI5, mxI6] = ind2sub(size(meas1), mxInd2);
fprintf('\nOverall Max:             %f',mxVal2)
fprintf('\n Overall Max Location: (%i, %i, %i, %i, %i, %i)', ...
        mxI1, mxI2, mxI3, mxI4, mxI5, mxI6);
% Take the max in the frame with the min
meas3 = meas1(:,:, I3, I4, I5, I6); %Frame with min
mxVal = max(meas3(:)); 
fprintf('\nMax in min Frame:        %f', mxVal)
fprintf('\nMax/Min, min Frame-dark: %f', mxVal/mnVal)
%fprintf('\nMax/Min ratio-dark:      %f\n', (mxVal-0.005)/(mnVal-0.005))
if ~(isPMNorm || isMMFNorm)
    fprintf('\nWARN: redPM not used; fib tip power provided by user');
end
eta_s = mnVal/pmRead1;
eta_p = mxVal/pmRead1;
fprintf('\nPower on Fiber tip:      %f', pmRead1);
fprintf('\nMin, PM normed (eta_s):  %f', eta_s);
fprintf('\nMax, PM normed (eta_p):  %f\n', eta_p);

if ~sum(strfind(expNm, 'PSF'))
    % Assume we are analyzing a donut
    if ~(isPMNorm || isMMFNorm)
        fprintf('\n      WARN: redPM not used; fib tip power provided by user');
    end
    fprintf('\n      - Rel. Int. Time = (%f ^2/%f) - bias corr.------> = %0.2f\n', eta_p, eta_s, eta_p*eta_p/eta_s)
    fprintf('      - Ratio in Frame = (%f/%f) - bias corr.------> = %0.2f\n', mxVal, mnVal, eta_p/eta_s)
    fprintf('      - Eta_s = (%f/%f) w/ approx vals for thrpt -> = %f\n', mnVal, pmRead1, eta_s)
    fprintf('      - Eta_p = (%f/%f) w/ approx vals for thrpt -> = %f\n', mxVal, pmRead1, eta_p)
    
    %Automated DataNotes
    fileID = fopen(datFl, 'a');
    fprintf(fileID, '\n         - Rel. Int. Time = (%f ^2/%f) - bias corr.------> = %0.2f\n', eta_p, eta_s, eta_p*eta_p/eta_s);
    fprintf(fileID, '         - Ratio in Frame = (%f/%f) - bias corr.------> = %0.2f\n', mxVal, mnVal, eta_p/eta_s);
    fprintf(fileID, '         - Eta_s = (%f/%f) w/ approx vals for thrpt -> = %f\n', mnVal, pmRead1, eta_s);
    fprintf(fileID, '         - Eta_p = (%f/%f) w/ approx vals for thrpt -> = %f\n\n', mxVal, pmRead1, eta_p);
    fclose(fileID);
else
    % Analysing the PSF
    fprintf(['\n Max Location: ', ...
        '\n  X index = %3i, X   = %6.2f     V', ...
        '\n  Y index = %3i, Y   = %6.2f     V', ...
        '\n  Foc ind = %3i, Foc =   %0.5f  mm', ...
        '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
        '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
        mxI1, distX(mxI1), mxI2, distY(mxI2), ...
        mxI4, distZ(mxI4), mxI5, distVX(mxI5), mxI6, distVY(mxI6));
    if ~(isPMNorm || isMMFNorm)
        fprintf('\n      WARN: redPM not used; fib tip power provided by user');
    end
    fprintf('\n      - Max:                        %f', mxVal2)
    fprintf('\n      - Power on Fiber tip:         %f',pmRead1)
    fprintf('\n      - Max, PM normed (coupling):  %f\n',mxVal2/pmRead1)
    % Automated DataNotes
    fileID = fopen(datFl, 'a');
    fprintf(fileID, '\r\n         - Max:                        %f', mxVal2);
    fprintf(fileID, '\r\n         - Power on Fiber tip:         %f',pmRead1);
    fprintf(fileID, '\r\n         - Max, PM normed (coupling):  %f\r\n\r\n',mxVal2/pmRead1);
    fclose(fileID);
    % Display max in each frame of focus
    tmp1 = nan(size(measScl,4),1);
    for iii = 1:size(measScl, 4)
        tmp1(iii) = max(max(measScl(:,:,:,iii)));
    end
    fprintf('\n Coupling at all focci:\n')
    disp(tmp1')
end
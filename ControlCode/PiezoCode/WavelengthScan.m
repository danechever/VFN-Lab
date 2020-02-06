% WavelengthScan conducts fiber scans through bandwidths and wavelengths to
% record how the null and coupling are affected by the bandwidth and
% wavelength of light

% When setting wavelengths to scan, the script takes the min and max values
% indicated and creates an even number of points between them, including
% the minimum and maximum values
% For example, 
% WLPoints = 6 and WLMinMax = [400,600] would scan
% 400,440,480,520,560,600nm

% When setting bandwidths, the script scans evenly spaced bandwidths
% between 0 and MaxBW not including zero, with the number of bandwidths
% indicated in ScanPoi 
% for example,
% MaxBW = 60 and ScanPoi = 10 would scan 6,12,18,24,30,36,42,48,54,60 nm
% bandwidth for every wavelength center indicated previously

% When there are multiple bandwidths and wavelenths to scan, the script
% will scan every combination of the two, so if there are 4 wavelengths and
% 4 bandwidths being scanned, there will be 16 measurements - be careful
% with this because some of the analysis code cannot handle this and I have
% never tested this feature

% This script is also capable of changing the focus axis of the fiber to
% account for defocus in a linear fashion

% NOTE: Vortex must be in position prior to running script

clear all; close all;
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode'));

global s FMTO_scale

%% Parameters to change
% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\020520_FLE6';

% Experiment name for figures
expNm = '5WLS_Var755805NullCheck1';
message = 'Check null depth at various wavelengths with current system params.';
% Figures being created
etaSPlot = true;
etaPPlot = true;

% Wavelength centers being scanned
WLPoints = 5;
% Minimum WL and Maximum WL being scanned [Min,Max]
WLMinMax = [755,805];

%Highest Bandwidth
MaxBW = 3;
%Number of points being scanned
ScanPoi = 1;
%Varia power values for normalization (these should be the 'raw' values, so no femto responsivity factor)
VarPow = [ 0.08996; 0.25004; 0.40561; 0.16661; 0.15985];
%For 790nm VarPow = [23.60209 42.40228 55.91581 67.03991 77.09343 84.71259 93.21918 98.27705 102.61063 106.61146];

%-- General Varia Settings
varPWR = 80;
varEMS = true;

%-- Fiber Scan settings
xCenter = 76.00;
yCenter = 87.00;
zCenter = 1.73700;
backlash  = 10; 
xPoints = 30;
yPoints = xPoints;
VXcenter = 12.011700; 
VYcenter = 12.860000;
Vbacklash = .005;
% NOTE: code assumes null is min in image so refStep and centering should
  % be s.t. this is true.
refStep = 5;

%-- Focus settings %1.742 for 750nm %1.702 for 830
isAutoFocus = false && length(WLCenters)>1;   %if true, changes focus during wavelength scan in linear fashion
focusRange = [1.737 1.737]; %Optimal focus goes as focusRange(1) with WLMinMax(1)

%isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
Nread   = 100;      % power samples at a given locaiton
Nrate   = 1000;     % rate of scan [samples/second]

%% DataNotes automation pt1
% Write the filename
datName = strcat('\DataNotes.txt');
datFl = strcat(svFld, datName);
fileID = fopen(datFl, 'a');
fprintf(fileID, '- %s\r\n', expNm);

% General comment
 fprintf(fileID, '     %s\r\n', message);
fprintf(fileID, '    * Acquire the new rough data for Spirit of Lyot.\r\n');

% Scan Params section
fprintf(fileID, '    Scan params:\r\n');  
trash = [VXcenter; VYcenter];
fprintf(fileID, '        Vortex pos:    x = %9.6fmm    y = %9.6fmm\r\n', trash);
fprintf(fileID, '        Vbacklash:     %5.3fmm\r\n', Vbacklash); 
trash = [varPWR, mean(WLMinMax), MaxBW/ScanPoi, MaxBW];
fprintf(fileID, '        Laser power:   Varia: %d pct power, %6.4fnm center, %3.1fnm-%dnm bandwidth\r\n', trash);
fprintf(fileID, '        Optical Filt:  600nm Longpass \r\n');
fprintf(fileID, '        Vortex Mask:   JPL Poly 550nm Charge 2\r\n');
trash = [xCenter; yCenter; zCenter];
fprintf(fileID, '        scan center    (%5.2f V,%6.2f V, %8.6f mm)\r\n', trash);
trash = [xPoints; 1; refStep; refStep; xPoints; 1];
fprintf(fileID, '        Scan params:   xpoi = %d, zpoi = %d, refStep = %d, stepsize = %d/(%d/10); zstep = %4.2f mm\r\n', trash);
fprintf(fileID, '        isAutoScale:   %s\r\n', mat2str(true));
trash = [FMTO_scale; Nread; Nrate];
fprintf(fileID, '        Femto:         FMTO_scale = %d; Nread = %d; Nrate = %d\r\n', trash);
% trash = [pmNormReport; pmCalWvl];
% fprintf(fileID, '        redPM norm:    isPMNorm = %s; pmNormReport = %5.2f; pmCalWavelength = %dnm\r\n\r\n', mat2str(isPMNorm), trash);

fclose(fileID);

%~~ END DataNotes pt1 automation
%% Initialize Varia
% Connect and instantiate NKT
VFN_setUpNKT;

%-- Set varia to desired initial setting
% Set power
VFN_NKT_setPowerLevel(NKT, varPWR);
% Set emission
VFN_NKT_setEmission(NKT, varEMS);

%% Initialize the Piezos
%Create serial port and define communication properties
    %COM# changes between computers; must check in Device Manager
MDT     = serial('COM4', 'BaudRate', 115200, 'InputBufferSize', 1500, ...
    'Terminator', {'CR' 'LF/CR'});
fopen(MDT);           %Connect port to device
%% Zaber Setup

zabbacklash= 0.005;           %Backlash removal step size for zabers

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

VFN_Zab_move(fibZ, zCenter-zabbacklash);    
VFN_Zab_move(fibZ, zCenter);
    
%% Femto setup
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

VFN_FMTO_LUCI_setGain(FMTO_scale);

fprintf('Current Gain setting: %i\n', FMTO_scale)
%% Main loop
%Estimate on scan time
fprintf('Estimated Time:    %6.1f minutes\n', WLPoints*ScanPoi*xPoints*xPoints*.5/60 + 1*ScanPoi*WLPoints);

% Wavelength centers being scanned
WLCenters = linspace(WLMinMax(1),WLMinMax(2), WLPoints);

% Define vector of bandwidths to sample
BWs = linspace(round(MaxBW/ScanPoi), MaxBW, ScanPoi);

if isAutoFocus
    % Define vector of foci
    foci = linspace(focusRange(1), focusRange(2), WLPoints);
end

% --Define raw data vectors
% Values of Nulls
dat = nan(length(WLCenters),length(BWs),1);
% Values of Peaks
datmax = nan(length(WLCenters),length(BWs),1);
% All data
dat2D = nan(length(WLCenters),length(BWs),xPoints + 1,yPoints + 1);
rawdat = nan(WLPoints, ScanPoi, xPoints + 1, yPoints + 1, Nread);
scls = nan(WLPoints, ScanPoi, xPoints + 1, yPoints + 1, 3);

% Define raster scan properties
stepSize = refStep/(xPoints/10);
distX = xCenter-stepSize*xPoints/2:stepSize:xCenter+stepSize*xPoints/2;
distY = yCenter-stepSize*yPoints/2:stepSize:yCenter+stepSize*yPoints/2;
dist = [distX', distY'];

% Loop through all wavelengths
for j = 1:length(WLCenters)
    if isAutoFocus
        % Change focus
        VFN_Zab_move(fibZ, foci(j)-zabbacklash);    
        VFN_Zab_move(fibZ, foci(j));
    end
    disp(['Wavelength set to ' num2str(WLCenters(j))]);
    
    % Loop through all bandwidths
    for i = 1:length(BWs)
        %-- Set Varia
        % Change bandwidth
        VFN_NKT_setWvlRange(NKT, WLCenters(j)-BWs(i)/2,WLCenters(j)+BWs(i)/2);
        %pause for 60 seconds to let Varia stabilize
        pause(60)
        disp(['Bandwidth set to ' num2str(BWs(i))]);

        %--Scan fiber
        % norm at current varia settings
        %pmNormReport = VarPow(j,i);   
        % Use courseRaster to do scan
        [measScl, rawdat(j,i,:,:,:), scls(j,i,:,:,:)] = CourseRaster_FullScan(MDT, dist);

        % Assume null is min in image
        dat(j,i) = min(measScl(:));
        datmax(j,i) = max(measScl(:));

        % Save 2D raster result matrix
        dat2D(j,i,:,:) = measScl;
    end
end
if size(VarPow) == size(dat2D(:,:,1,1))
    for j = 1:length(WLCenters)
        eta_s(j,:) = dat(j,:) ./ (VarPow*VFN_getFmtoResponsivity(WLCenters(j)));
        eta_p(j,:) = datmax(j,:) ./ (VarPow*VFN_getFmtoResponsivity(WLCenters(j)));
    end
end

%% Clean-Up devices
%-- Clean up NKT
VFN_cleanUpNKT;

% Set Femto back to 10^6 
FMTO_scale = 6;
VFN_FMTO_LUCI_setGain(FMTO_scale);
fprintf('\n Femto Gain set to %i\n',FMTO_scale);

%Set all three axes to starting position
DE2_MDTVol(MDT, xCenter-backlash, 'x', 0); %xC holds current xVoltage
DE2_MDTVol(MDT, xCenter, 'x', 0); %xC holds current xVoltage
DE2_MDTVol(MDT, yCenter-backlash, 'y', 0); %yC holds current yVoltage
DE2_MDTVol(MDT, yCenter, 'y', 0); %yC holds current yVoltage
% DE2_MDTVol(MDT, Zcenter-backlash, 'z', 0); %zC holds current zVoltage
% DE2_MDTVol(MDT, Zcenter, 'z', 0); %zC holds current zVoltage

%if isAutoFocus %Change by BC: isAutoFocus isn't checked to init zabs
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
%end
fclose(MDT);
delete(MDT);
clear MDT;
%% Display Data
for b = 1:length(WLCenters)
    for a = 1:length(BWs)
        figure();
        colormap('gray');
        %Not sure why transpose and "axis xy" are needed for correct orient.
        imagesc((distX-xCenter),(distY-yCenter),transpose(squeeze(dat2D(b,a,:,:))));
        titStr = sprintf('PWR Thru Fib');
        title(titStr)
        axis xy;axis image%axis square;
        cbr = colorbar;
        ylabel(cbr, 'Power [V]')
        set(gca,'TickDir','out')
        xlabel('Fiber X position [V]')
        ylabel('Fiber Y position [V]')
        tag = sprintf('%ix%i',xPoints,yPoints);
        tag = sprintf([tag 'Wvl%i_BW%i'],WLCenters(b), BWs(a));

        saveas(gcf, [svFld filesep expNm '_' tag '_CoupMap.png'])
    end
end
% etaS plots
if etaSPlot == true
    if length(BWs)>1
        for j = 1:length(WLCenters)
            figure();
            scatter(BWs,dat(j,:))
            titStr = sprintf('Power Vs Bandwidth');
            title(titStr)
            xlabel(['Varia Bandwidth about ' num2str(WLCenters(j)) 'nm'])
            ylabel('Power')
            saveas(gcf, [svFld filesep expNm '_BWGraph_EtaS_CoupMap.png'])
            if size(VarPow) == size(dat2D(:,:,1,1))
                figure();
                scatter(BWs,eta_s(j,:))
                titStr = sprintf('Eta_s Vs Bandwidth');
                title(titStr)
                xlabel(['Varia Bandwidth about ' num2str(WLCenters(j)) 'nm'])
                ylabel('\eta_{s}')
                saveas(gcf, [svFld filesep expNm '_BWGraph_EtaS_CoupMap_Norm.png'])
            end
        end
    end

    if length(WLCenters)>1
        for i = 1:length(BWs)
            figure();
            scatter(WLCenters,dat(:,i))
            titStr = sprintf('Power Vs Wavelength');
            title(titStr)
            xlabel(['Varia Wavelength at ' num2str(BWs(i)) 'nm BW'])
            ylabel('Power')
            saveas(gcf, [svFld filesep expNm '_WvlGraph_EtaS_CoupMap.png'])
            if size(VarPow) == size(dat2D(:,:,1,1))
                    figure();
                    scatter(WLCenters,eta_s(:,i))
                    titStr = sprintf('Eta_s Vs Wavelength');
                    title(titStr)
                    xlabel(['Varia Wavelength at ' num2str(BWs(i)) 'nm BW'])
                    ylabel('\eta_{s}')
                    saveas(gcf, [svFld filesep expNm '_WvlGraph_EtaS_CoupMap_Norm.png'])
            end
        end
    end
end
% etaP plots
if etaPPlot == true
    if length(BWs)>1
        for j = 1:length(WLCenters)
            figure();
            scatter(BWs,dat(j,:))
            titStr = sprintf('Power Vs Bandwidth');
            title(titStr)
            xlabel(['Varia Bandwidth about ' num2str(WLCenters(j)) 'nm'])
            ylabel('Power')
            saveas(gcf, [svFld filesep expNm '_BWGraph_EtaP_CoupMap.png'])
            if size(VarPow) == size(dat2D(:,:,1,1))
                figure();
                scatter(BWs,eta_p(j,:))
                titStr = sprintf('Eta_p Vs Bandwidth');
                title(titStr, 'interpreter', 'none')
                xlabel(['Varia Bandwidth about ' num2str(WLCenters(j)) 'nm'])
                ylabel('\eta_{p}')
                saveas(gcf, [svFld filesep expNm '_BWGraph_EtaP_CoupMap_Norm.png'])
            end
        end
    end

    if length(WLCenters)>1
        for i = 1:length(BWs)
            figure();
            scatter(WLCenters,dat(:,i))
            titStr = sprintf('Power Vs Wavelength');
            title(titStr)
            xlabel(['Varia Wavelength at ' num2str(BWs(i)) 'nm BW'])
            ylabel('Power')
            saveas(gcf, [svFld filesep expNm '_WvlGraph_EtaP_CoupMap.png'])
            if size(VarPow) == size(dat2D(:,:,1,1))
                    figure();
                    scatter(WLCenters,eta_p(:,i))
                    titStr = sprintf('Eta_p Vs Wavelength');
                    title(titStr, 'interpreter', 'none')
                    xlabel(['Varia Wavelength at ' num2str(BWs(i)) 'nm BW'])
                    ylabel('\eta_{p}')
                    saveas(gcf, [svFld filesep expNm '_WvlGraph_EtaP_CoupMap_Norm.png'])
            end
        end
    end
end

%% Report results
disp('Null points')
if(length(BWs)>1)
    for j = 1:length(WLCenters)
        fprintf('\n        WvlCenter: %3dnm\n',WLCenters(j))
        for i = 1:length(BWs)
            if size(VarPow) == size(dat2D(:,:,1,1)) 
                fprintf('        BW: %2dnm, Null = %f, Coupling = %f\n',BWs(i), dat(j,i)/VarPow(j,i),datmax(j,i)/VarPow(j,i))
            else
                fprintf('        BW: %2dnm, Null (Unnorm) = %f, Coupling (Unnorm) = %f\n',BWs(i), dat(j,i), datmax(j,i))
            end
        end
    end
else
    fprintf('\n        BW: %2dnm\n',BWs(1))
    for j = 1:length(WLCenters)
        fprintf('        WvlCenter: %3dnm, Null = %f, Coupling = %f\n',WLCenters(j), dat(j,i)/VarPow(j,i),datmax(j,i)/VarPow(j,i))
    end
end

% %-- Use user-provided norm value if desired
% if ~isPMNorm
%     pmRead1 = pmNormReport;
% end
% 
% % Take average of all the Nread
% meas1 = measScl;
% % Crop the matrix to middle (cropVal-2)/cropVal region to be in donut null for min
% cropVal = 7;
% xmin = max(floor(length(distX)/cropVal),1);
% ymin = max(floor(length(distY)/cropVal),1);
% xmax = (floor((cropVal-1)*length(distX)/cropVal+1));
% ymax = (floor((cropVal-1)*length(distY)/cropVal+1));
% meas2 = meas1(xmin:xmax,ymin:ymax,:,:,:,:);
% 
% % Take the min of the cropped region along all dimensions 
% %mnVal= min(min(min(meas1(Xpoints/5:4*Xpoints/5+1,Ypoints/5:4*Ypoints/5+1,:))));
% [mnVal, mnInd] = min(meas2(:)); %get min and linear index
% [I1, I2, I3, I4, I5, I6] = ind2sub(size(meas2), mnInd); %linear ind to subind
% fprintf('\nMin in center region:    %f',mnVal)
% fprintf(['\n Min Location: ', ...
%          '\n  X index = %3i, X   = %6.2f     V', ...
%          '\n  Y index = %3i, Y   = %6.2f     V', ...
%          '\n  Foc ind = %3i, Foc =   %0.5f  mm', ...
%          '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
%          '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
%          I1+xmin-1, distX(I1+xmin-1), I2+ymin-1, distY(I2+ymin-1), ...
%          I4, distZ(I4), I5, distVX(I5), I6, distVY(I6));
% 
% 
% %Take overall max
% %mxVal = max(max(max(meas1)));
% %ind = floor(find(meas1>= mxVal)/((Xpoints+1)*(Ypoints+1)));
% %fprintf('\nFrame for max:        %d',ind)
% [mxVal2, mxInd2] = max(meas1(:));
% [mxI1, mxI2, mxI3, mxI4, mxI5, mxI6] = ind2sub(size(meas1), mxInd2);
% fprintf('\nOverall Max:             %f',mxVal2)
% fprintf('\n Overall Max Location: (%i, %i, %i, %i, %i, %i)', ...
%         mxI1, mxI2, mxI3, mxI4, mxI5, mxI6);
% % Take the max in the frame with the min
% meas3 = meas1(:,:, I3, I4, I5, I6); %Frame with min
% mxVal = max(meas3(:)); 
% fprintf('\nMax in min Frame:        %f', mxVal)
% fprintf('\nMax/Min, min Frame-dark: %f', mxVal/mnVal)
% %fprintf('\nMax/Min ratio-dark:      %f\n', (mxVal-0.005)/(mnVal-0.005))
% if ~isPMNorm
%     fprintf('\nWARN: redPM not used; fib tip power provided by user');
% end
% fprintf('\nPower on Fiber tip:      %f', pmRead1);
% fprintf('\nMin, PM normed (eta_s):  %f', mnVal/pmRead1);
% fprintf('\nMax, PM normed (eta_p):  %f\n', mxVal/pmRead1);
% 
% if ~sum(strfind(expNm, 'PSF'))
%     % Assume we are analyzing a donut
%     if ~isPMNorm
%         fprintf('\n      WARN: redPM not used; fib tip power provided by user');
%     end
%     fprintf('\n      - Ratio in Frame = (%f/%f) - bias corr.------> = %0.2f\n', mxVal, mnVal, mxVal/mnVal)
%     fprintf('      - Eta_s = (%f/%f) w/ approx vals for thrpt -> = %f\n', mnVal, pmRead1, mnVal/pmRead1)
%     fprintf('      - Eta_p = (%f/%f) w/ approx vals for thrpt -> = %f\n', mxVal, pmRead1, mxVal/pmRead1)
% else
%     % Analysing the PSF
%     fprintf(['\n Max Location: ', ...
%         '\n  X index = %3i, X   = %6.2f     V', ...
%         '\n  Y index = %3i, Y   = %6.2f     V', ...
%         '\n  Foc ind = %3i, Foc =   %0.5f  mm', ...
%         '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
%         '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
%         mxI1, distX(mxI1), mxI2, distY(mxI2), ...
%         mxI4, distZ(mxI4), mxI5, distVX(mxI5), mxI6, distVY(mxI6));
%     if ~isPMNorm
%         fprintf('\n      WARN: redPM not used; fib tip power provided by user');
%     end
%     fprintf('\n      - Max:                        %f', mxVal2)
%     fprintf('\n      - Power on Fiber tip:         %f',pmRead1)
%     fprintf('\n      - Max, PM normed (coupling):  %f\n',mxVal2/pmRead1)
%     % Display max in each frame of focus
%     tmp1 = nan(size(measScl,4),1);
%     for iii = 1:size(measScl, 4)
%         tmp1(iii) = max(max(measScl(:,:,:,iii)));
%     end
%     fprintf('\n Coupling at all focci:\n')
%     disp(tmp1')
% end
%% fits save
import matlab.io.*
tag = sprintf('fullCube');
% Cube names
datCubeNm   = [svFld filesep expNm '_' tag '.fits'];
GnCubeNm    = [svFld filesep expNm '_' tag '_GAINS.fits'];
SRCubeNm    = [svFld filesep expNm '_' tag '_SEMIREDU.fits'];
NormCubeNm    = [svFld filesep expNm '_' tag '_NORMREDU.fits'];

% Define permute order:
permOrd = [4,3,5,1,2];

%-- SAVE RAW DATA CUBE
fitmap= fits.createFile(datCubeNm);
fits.createImg(fitmap,'double',size(permute(rawdat,permOrd)));
% %header data
fits.writeKey(fitmap, 'Xcenter ', xCenter);
fits.writeKey(fitmap, 'Xpoints ', xPoints);
fits.writeKey(fitmap, 'Ycenter ', yCenter);
fits.writeKey(fitmap, 'Ypoints ', yPoints);
fits.writeKey(fitmap, 'XYsteps ', stepSize);
%fits.writeKey(fitmap, 'Zsteps ',  ZStepSize);
fits.writeKey(fitmap, 'WvlMin ', WLMinMax(1));
fits.writeKey(fitmap, 'WvlMax ', WLMinMax(2));
fits.writeKey(fitmap, 'WvlPts', WLPoints);
fits.writeKey(fitmap, 'BWMin ', MaxBW/ScanPoi);
fits.writeKey(fitmap, 'BWMax ', MaxBW);
fits.writeKey(fitmap, 'BWPts', ScanPoi);
fits.writeKey(fitmap, 'Nread ',  Nread);
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', 'Nread');
fits.writeKey(fitmap, 'NAX4', 'Wavelength');
fits.writeKey(fitmap, 'NAX5', 'Bandwidth');
%Zaber status
% fits.writeKey(fitmap, 'isZab ',  logical(isZab));
% fits.writeKey(fitmap, 'VXcenter',  VXcenter);
% fits.writeKey(fitmap, 'VYcenter',  VYcenter);
% fits.writeKey(fitmap, 'isVortSc', logical(isVortScan));
% fits.writeKey(fitmap, 'VXpoints', VXpoints);
% fits.writeKey(fitmap, 'VYpoints', VYpoints);
% fits.writeKey(fitmap, 'vStepRng', vStepRange);
% Power meter level
fits.writeKey(fitmap, 'isAutScl', logical(true));
ind = strfind(GnCubeNm ,'\');
fits.writeKey(fitmap, 'DatCube ', GnCubeNm(ind(end-1)+1:end));
ind = strfind(SRCubeNm ,'\');
fits.writeKey(fitmap, 'SRCube ', SRCubeNm(ind(end-1)+1:end));
ind = strfind(NormCubeNm ,'\');
fits.writeKey(fitmap, 'NormCube ', NormCubeNm(ind(end-1)+1:end));
%fits.writeKey(fitmap, 'NormFlag', logical(isPMNorm));
%fits.writeKey(fitmap, 'NormMean', mean(pmRead));
% Saving Normalization data
for j = 1:length(WLCenters)
    for i = 1:ScanPoi
        fits.writeKey(fitmap, ['RAWNRM' sprintf('%02d', (j-1)*ScanPoi+i)], VarPow(j,i));
    end
end
for j = 1:length(WLCenters)
    for i = 1:ScanPoi
        fits.writeKey(fitmap, ['NRMVAL' sprintf('%02d', (j-1)*ScanPoi+i)], VarPow(j,i)*VFN_getFmtoResponsivity(WLCenters(j)));
    end
end

fits.writeKey(fitmap, 'CNTRLCD', 'WavelengthScan');
% 
fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(rawdat,permOrd));
fits.closeFile(fitmap);

%-- SAVE GAIN SETTINGS
fitmap= fits.createFile(GnCubeNm);
fits.createImg(fitmap,'double',size(permute(scls, permOrd)));
%header data
ind = strfind(datCubeNm,'\');
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', '1=FMTO_Scale;2=Gain;3=Bias');
fits.writeKey(fitmap, 'NAX4', 'Wavelengths');
fits.writeKey(fitmap, 'NAX5', 'Bandwidth');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(scls,permOrd));
fits.closeFile(fitmap);

%-- SAVE SEMI REDUCED DATA
permOrd = [4,3,1,2];
fitmap= fits.createFile(SRCubeNm);
fits.createImg(fitmap,'double',size(permute(dat2D, permOrd)));
%header data
ind = strfind(datCubeNm,'\');
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
fits.writeKey(fitmap, 'NAX3', 'Wavelengths');
fits.writeKey(fitmap, 'NAX4', 'Bandwidth');
%fits.writeKey(fitmap, 'NAX4', 'Focus');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(dat2D,permOrd));
fits.closeFile(fitmap);
 
if (size(VarPow) == size(dat2D(:,:,1,1)))
    dat2Dnorm = nan(size(dat2D));
    for j = 1:WLPoints
        for i = 1:ScanPoi
            dat2Dnorm(j,i,:,:) = dat2D(j,i,:,:) ./ VarPow(j,i);
        end
    end
    
    
    %-- SAVE SEMI REDUCED NORM'D DATA
    fitmap= fits.createFile(NormCubeNm);
    fits.createImg(fitmap,'double',size(permute(dat2Dnorm,permOrd)));
    %header data
    ind = strfind(datCubeNm,'\');
    fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
    fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
    fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
    fits.writeKey(fitmap, 'NAX3', 'Wavelengths');
    fits.writeKey(fitmap, 'NAX4', 'Bandwidth');
    %fits.writeKey(fitmap, 'NAX4', 'Focus');
    %fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
    %fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');
    %fits.writeKey(fitmap, 'NormFlag', logical(isPMNorm));
    %fits.writeKey(fitmap, 'NormMean', mean(pmRead));
    %fits.writeKey(fitmap, 'NormSTD', std(pmRead));

    fits.setCompressionType(fitmap,'NOCOMPRESS');
    fits.writeImg(fitmap,permute(dat2Dnorm,permOrd));
    fits.closeFile(fitmap);
end


clear all; close all;
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\VFN-Lab\ControlCode'));

global s FMTO_scale

%% Parameters to change
% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\071019_VAR1';

% Experiment name for figures
expNm = '13BWS_Var640BWScan2';

% Wavelength centers being scanned
WLPoints = 1;
% Minimum WL and Maximum WL being scanned
WLMinMax = [635,800];

%Highest Bandwidth
MaxBW = 60;
%Number of points being scanned
ScanPoi = 10;
%Varia power values for normalization
VarPow = [.108 .108 .108 .108 .108];

%-- General Varia Settings
varPWR = 80;
varEMS = true;

%-- Fiber Scan settings
xCenter = 118.15;
yCenter = 81.72;
backlash  = 10; 
xPoints = 50;
yPoints = xPoints;
% NOTE: code assumes null is min in image so refStep and centering should
  % be s.t. this is true.
refStep = 5;

isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
Nread   = 100;      % power samples at a given locaiton
Nrate   = 1000;     % rate of scan [samples/second]

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

%% Femto setup
% Define the gain to apply for scaling. Uses same gain for all: gnFact^(FMTO_scale-6)
gnFact  = 9.97;

% Setup Femto 
VFN_setUpFMTO;  

s = VFN_FMTO_setGain(s, FMTO_scale);

fprintf('Current Gain setting: %i\n', FMTO_scale)
%% Main loop
%Estimate on scan time
fprintf('Estimated Time:    %6.1f minutes\n', WLPoints*ScanPoi*xPoints*xPoints*.5/60 + .5*ScanPoi);

% Wavelength centers being scanned
WLCenters = linspace(WLMinMax(1),WLMinMax(2), WLPoints);

%Define vector of bandwidths to sample
BWs = linspace(round(MaxBW/ScanPoi), MaxBW, ScanPoi);

% --Define raw data vectors
% Values of Nulls
dat = nan(length(WLCenters),length(BWs),1);
% All data
dat2D = nan(length(WLCenters),length(BWs),xPoints + 1,yPoints + 1);

% Define raster scan properties
stepSize = refStep/(xPoints/10);
distX = xCenter-stepSize*xPoints/2:stepSize:xCenter+stepSize*xPoints/2;
distY = yCenter-stepSize*yPoints/2:stepSize:yCenter+stepSize*yPoints/2;
dist = [distX', distY'];

%Loop through all wavelengths in Bandwidth
for j = 1:length(WLCenters)
    for i = 1:length(BWs)
        %-- Set Varia
        % Change bandwidth
        VFN_NKT_setWvlRange(NKT, WLCenters(j)-BWs(i)/2,WLCenters(j)+BWs(i)/2);
        %pause for 60 seconds to let Varia stabilize
        pause(60)
        disp(['Bandwidth set to ' num2str(BWs(i))]);

        %--Scan fiber
        % norm at current varia settings
        pmNormReport = VarPow(j,i);   
        % Use courseRaster to do scan
        measScl = CourseRaster_FullScan(MDT, dist);

        % Assume null is min in image
        dat(j,i) = min(measScl(:));

        % Save 2D raster result matrix
        dat2D(j,i,:,:) = measScl;
    end
end

eta_s = dat./VarPow';

%% Clean-Up devices
%-- Clean up NKT
VFN_cleanUpNKT;

% Set Femto back to 10^6 
FMTO_scale = 6;
s = VFN_FMTO_setGain(s, FMTO_scale);
fprintf('\n Femto Gain set to %i',FMTO_scale);

%Set all three axes to starting position
DE2_MDTVol(MDT, xCenter-backlash, 'x', 0); %xC holds current xVoltage
DE2_MDTVol(MDT, xCenter, 'x', 0); %xC holds current xVoltage
DE2_MDTVol(MDT, yCenter-backlash, 'y', 0); %yC holds current yVoltage
DE2_MDTVol(MDT, yCenter, 'y', 0); %yC holds current yVoltage
% DE2_MDTVol(MDT, Zcenter-backlash, 'z', 0); %zC holds current zVoltage
% DE2_MDTVol(MDT, Zcenter, 'z', 0); %zC holds current zVoltage

fclose(MDT);
delete(MDT);
clear MDT;
%% Display Data
for b = 1:length(WLCenters)
    for a = 1:length(BWs)
        figure();
        colormap('gray');
        %Not sure why transpose and "axis xy" are needed for correct orient.
        imagesc((distX-xCenter),(distY-yCenter),transpose(squeeze(dat2D(a,:,:))));
        titStr = sprintf('PWR Thru Fib');
        title(titStr)
        axis xy;axis image%axis square;
        cbr = colorbar;
        ylabel(cbr, 'Power [V]')
        set(gca,'TickDir','out')
        xlabel('Fiber X position [V]')
        ylabel('Fiber Y position [V]')
        tag = sprintf('%ix%i',xPoints,yPoints);
        tag = sprintf([tag '_BW%i'],BWs(a));

        saveas(gcf, [svFld filesep expNm '_' tag '_CoupMap.png'])
    end
end
if length(BWs)>1
    for j = 1:length(WLCenters)
        figure();
        scatter(BWs,dat(j,:))
        titStr = sprintf('Power Vs Bandwidth');
        title(titStr)
        xlabel(['Varia Bandwidth about ' num2str(WLCenters(j)) 'nm'])
        ylabel('Power')
        saveas(gcf, [svFld filesep expNm '_BWGraph_CoupMap.png'])

        figure();
        scatter(BWs,eta_s(j,:))
        titStr = sprintf('Eta_s Vs Bandwidth');
        title(titStr)
        xlabel(['Varia Bandwidth about ' num2str(WLCenters(j)) 'nm'])
        ylabel('Eta s')
        saveas(gcf, [svFld filesep expNm '_BWGraph_CoupMap_Norm.png'])
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
        saveas(gcf, [svFld filesep expNm '_WvlGraph_CoupMap.png'])

        figure();
        scatter(WLCenters,eta_s(:,i))
        titStr = sprintf('Eta_s Vs Wavelength');
        title(titStr)
        xlabel(['Varia Wavelength at ' num2str(BWs(i)) 'nm BW'])
        ylabel('Eta s')
        saveas(gcf, [svFld filesep expNm '_WvlGraph_CoupMap_Norm.png'])
    end
end

%% Report results
disp('Null points')
if(length(BWs)>1)
    for j = 1:length(WLCenters)
        fprintf('\nWvlCenter: %3dnm\n',WLCenters(j))
        for i = 1:length(BWs)
            fprintf('BW: %2dnm, Null = %f\n',BWs(i), dat(j,i))
        end
    end
else
    fprintf('\nBW: %2dnm\n',BWs(1))
    for j = 1:length(WLCenters)
        fprintf('WvlCenter: %3dnm, Null = %f\n',WLCenters(j), dat(j,1))
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

% Define permute order:
permOrd = [3,2,1];

%-- SAVE DATA CUBE
fitmap= fits.createFile(datCubeNm);
fits.createImg(fitmap,'double',size(permute(dat2D,permOrd)));
% %header data
% fits.writeKey(fitmap, 'Xcenter ', xCenter);
% fits.writeKey(fitmap, 'Xpoints ', Xpoints);
% fits.writeKey(fitmap, 'Ycenter ', yCenter);
% fits.writeKey(fitmap, 'Ypoints ', Ypoints);
% fits.writeKey(fitmap, 'Zcenter ', Zcenter);
% %fits.writeKey(fitmap,'Zstart ',Zstart);
% fits.writeKey(fitmap, 'Zpoints ', Zpoints);
% fits.writeKey(fitmap, 'XYsteps ', StepSize);
% fits.writeKey(fitmap, 'Zsteps ',  ZStepSize);
% fits.writeKey(fitmap, 'Nread ',  Nread);
% fits.writeKey(fitmap, 'NAX1', 'XFiberCoord');
% fits.writeKey(fitmap, 'NAX2', 'YFiberCoord');
% fits.writeKey(fitmap, 'NAX3', 'Nread');
% fits.writeKey(fitmap, 'NAX4', 'Focus');
% fits.writeKey(fitmap, 'NAX5', 'XVortexCoord');
% fits.writeKey(fitmap, 'NAX6', 'YVortexCoord');
% %Zaber status
% fits.writeKey(fitmap, 'isZab ',  logical(isZab));
% fits.writeKey(fitmap, 'VXcenter',  VXcenter);
% fits.writeKey(fitmap, 'VYcenter',  VYcenter);
% fits.writeKey(fitmap, 'isVortSc', logical(isVortScan));
% fits.writeKey(fitmap, 'VXpoints', VXpoints);
% fits.writeKey(fitmap, 'VYpoints', VYpoints);
% fits.writeKey(fitmap, 'vStepRng', vStepRange);
% % Power meter level
% fits.writeKey(fitmap, 'isAutScl', logical(isAutoScale));
% ind = strfind(GnCubeNm ,'\');
% fits.writeKey(fitmap, 'DatCube ', GnCubeNm(ind(end-1)+1:end));
% ind = strfind(SRCubeNm ,'\');
% fits.writeKey(fitmap, 'SRCube ', SRCubeNm(ind(end-1)+1:end));
% ind = strfind(NormCubeNm ,'\');
% fits.writeKey(fitmap, 'NormCube ', NormCubeNm(ind(end-1)+1:end));
% fits.writeKey(fitmap, 'NormFlag', logical(isPMNorm));
% fits.writeKey(fitmap, 'NormMean', mean(pmRead));
% fits.writeKey(fitmap, 'NormSTD', std(pmRead));
% fits.writeKey(fitmap, 'CNTRLCD', 'PZScan_ZabFoc');
% 
fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,permute(dat2D,permOrd));
fits.closeFile(fitmap);

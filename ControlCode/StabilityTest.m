%{
    Code for commanding all the devices to a certain position and measuring
        the stability of the system at that point.
%}

%% add the necessary functions
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\ControlCode\'));

close all
clear all

%% General settings 

% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\102018_FTO5';

% Experiment name for figures
expNm = '16Sability_30MinRun_20sInterval';

% Number of frames/measurements to take
NMeas = 90;

% Interval between frames/measurements
Interv= 20;     %[s]

%~~ ZABER STUFF 
% Flag to control zaber motion (including centering)
isZab   = true;     % This controls the vortex zabers
isPMRd  = true;     % This controls the PM zaber and commands it to do a measure

% Vortex position properties [in mm]
VXcenter = 12.383445;      % Value when PM stage was added: 12.271445mm
VYcenter = 13.106083;  %!! must be >13.134083-0.8
Vbacklash= 0.005;
%~~ END ZABER STUFF 

%~~ PIEZO STUFF 
% Fiber position in Volts
Xcenter = 79.20;
Ycenter = 91.19;
Zcenter = 62.50;
backlash  = 10; % Kept backlash in to see if it did anything
%~~ END PIEZO STUFF 

%~~ POWER METER STUFF 
isAutoScale = true; % If false, the gain is held fixed. Else, gain
                        % is set automatically using FMTO_setAutoGain()
FMTO_scale = 6;     % Starting gain. the n in: 10^n for the gain setting
Nread   = 100;      % power samples at a given locaiton
Nrate   = 1000;     % rate of scan [samples/second]
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)
%~~ END POWER METER STUFF

%% Zaber Setup
if isZab || isPMRd
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
        pmX = ax93;
        clear ax93
    end
    
    % Move red PM out of beam
    VFN_Zab_move(pmX, 25);
    
end

%% MDT (piezo) Setup
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

%% Red Thorlabs PM setup
if isPMRd
    % Find a VISA-USB object.
    obj1 = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x1313::0x8078::P0015560::0::INSTR', 'Tag', '');
    
    % Create the VISA-USB object if it does not exist
    % otherwise use the object that was found.
    if isempty(obj1)
        obj1 = visa('NI', 'USB0::0x1313::0x8078::P0015560::0::INSTR');
    else
        fclose(obj1);
        obj1 = obj1(1);
    end
    
    % Connect to instrument object, obj1.
    fopen(obj1);
end

%% Perform scan
%Preallocate data matrices
meas=zeros(NMeas,Nread); 
measScl = meas;         %Matrix for semi-processed data
scales  = nan(NMeas,5);   %Matrix for scales, gains, biases, PMread, and time

if isZab
	%remove zaber X backlash
    VFN_Zab_move(vortX, VXcenter-Vbacklash);
    %move vortex to new X position
    fprintf('\nVortX pos: %f',VFN_Zab_move(vortX, VXcenter));
    %remove zaber y backlash
    VFN_Zab_move(vortY, VYcenter-Vbacklash);
    %move vortex to new Y position
    fprintf('  VortY pos: %f\n',VFN_Zab_move(vortY, VYcenter));
end
% Move Focus into place
DE2_MDTVol(MDT, Zcenter-backlash, 'z', 0);
DE2_MDTVol(MDT, Zcenter, 'z', 0);
% Move fiber X into place
DE2_MDTVol(MDT, Xcenter-backlash, 'x', 0);
DE2_MDTVol(MDT, Xcenter, 'x', 0);
% Move fiber Y into place
DE2_MDTVol(MDT, Ycenter-backlash, 'y', 0);
DE2_MDTVol(MDT, Ycenter, 'y', 0);
% Wait for femto power to settle
pause(Delay);
for i=1:NMeas
    %-- Iterate through measurements
    if mod(i,3) == 0
        fprintf(' Iter: %i \n',i)
    end

    % Mark when the measurement started
    tic;
    % Take a sample reading at current power
    read    = startForeground(s);
    if isAutoScale
        % Save old FMTO_scale for comparison after autoGain
        old_scale = FMTO_scale;
        % Save old read value for comparison after autoGain
        old_read = mean(read);
        % Modify gain accordingly
        [FMTO_scale, s] = VFN_FMTO_setAutoGain(s, old_read, FMTO_scale);
        % Check if re-read is needed
        if old_scale ~= FMTO_scale
            % FMTO_scale changed so the gain changed
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
            read    = startForeground(s);
            %ratio   = ratio * (old_read/mean(read));
        end
    end
    if find(read>10)
        warning('Power is too high')
    end
    meas(i,:) = read;
    scales(i,1) = FMTO_scale;
    scales(i,2) = gnFact^-(FMTO_scale-6);
    switch FMTO_scale
        case 5
            locBias = 0.0033887;
        case 6
            locBias = 0.0042960;
        case 7
            locBias = 0.0044690;
        case 8
            locBias = 0.0047430;
        case 9
            locBias = 0.0068927;
        otherwise
            warning('No bias forFMTO_scale = %i\n',FMTO_scale)
            locBias = nan;
    end
    scales(i,3) = locBias;
    if isPMRd
        if VFN_Zab_getPos(vortY) < 12.35
            error('Vortex Y-Axis is in the way')
        end
        % Move into beam if no collision
        VFN_Zab_move(pmX, 0);
        
        % Wait for PM value to settle
        pause(1);
        
        % Take 100 measurements of the power
        pmNread = 100;
        pmRead = nan(pmNread,1);
        for ii = 1:pmNread
            pmRead(ii)=str2num(query(obj1, 'measure:power?'));
        end
                
        % Move out of beam
        VFN_Zab_move(pmX, 25);
        scales(i,4) = mean(pmRead(:));
    end
    scales(i,5) = now;
    measScl(i,:) = read - locBias; %Subtract bias for semi-processed matrix
    
    % Wait for interval time to end
    while toc < Interv
        pause(0.25)
    end    
end

%% close the connections to all devices
%-- MDT
fclose(MDT);
delete(MDT);
clear MDT;

%-- Zabers
if isZab || isPMRd
    VFN_cleanUpZabers;
    if exist('vortX', 'var')
        clear vortX
    end
    if exist('vortY', 'var')
        clear vortY
    end
    if exist('pmX', 'var')
        clear pmX
    end
end

%-- Thorlabs PM
if isPMRd
    fclose(obj1);
end

%% Scale data accordingly to compensate for varied gain
%-- This is for plotting and analysis; the raw data is still stored as well

%** Bias is already subtracted in the main for-loop section

%-- Define scaled matrix; equal to mean of meas unless gain was varied
measScl = mean(meas,2).*scales(:,2);

%-- Define a separate matrix which normalizes the data by power on the fiber
% Account for lens losses
pmRead = scales(:,4)*0.98;
% Convert the calibration PM values from W to uW
pmRead = pmRead*10^6;
% Translate the calibration PM value from uW to V (at 10^6)
pmRead = pmRead*0.63;     % V/uW conversion is ~ 0.63
% Create the matrix and normalize by total power on fib tip
measScl2 = measScl./pmRead;

%% Display the result
%-- Figure of Femto data
figure();
plot(scales(:,5), measScl2);
hold on
scatter(scales(:,5),measScl2);
hold off
datetick('x', 'HH:MM:SS')
title(sprintf('Power Thru Fib: (%0.2f V, %0.2f V, %0.2f V) (%f , %f)', ...
                Xcenter,Ycenter,Zcenter, VXcenter, VYcenter));
set(gca,'TickDir','out')
xlabel('Time HH:MM:SS ')
ylabel('Normalized Power [V]')
saveas(gcf, [svFld filesep expNm '_NormedMeas.png'])

%-- Figure of PM Read data
figure();
%Not sure why transpose and "axis xy" are needed for correct orient.
plot(scales(:,5), pmRead);
hold on
scatter(scales(:,5),pmRead);
hold off
datetick('x', 'HH:MM:SS')
title(sprintf('Power At Fib Tip: (%0.2f V, %0.2f V, %0.2f V) (%f , %f)', ...
                Xcenter,Ycenter,Zcenter, VXcenter, VYcenter));
set(gca,'TickDir','out')
xlabel('Time HH:MM:SS ')
ylabel('Normalized Power from Red Thorlabs PM [V]')
saveas(gcf, [svFld filesep expNm '_PMValues.png'])

%% fits save
import matlab.io.*
% Cube names
datCubeNm   = [svFld filesep expNm '.fits'];
GnCubeNm    = [svFld filesep expNm '_SCALES.fits'];
SRCubeNm    = [svFld filesep expNm '_SEMIREDU.fits'];
NormCubeNm    = [svFld filesep expNm '_NORMREDU.fits'];

%-- SAVE DATA CUBE
fitmap= fits.createFile(datCubeNm);
fits.createImg(fitmap,'double',size(meas));
%header data
fits.writeKey(fitmap, 'Xcenter ', Xcenter);
fits.writeKey(fitmap, 'Ycenter ', Ycenter);
fits.writeKey(fitmap, 'Zcenter ', Zcenter);
fits.writeKey(fitmap, 'Nread ',  Nread);
fits.writeKey(fitmap, 'NAX1', 'Sample');
fits.writeKey(fitmap, 'NAX2', 'Nread');
%Zaber status
fits.writeKey(fitmap, 'isZab ',  logical(isZab));
fits.writeKey(fitmap, 'VXcenter',  VXcenter);
fits.writeKey(fitmap, 'VYcenter',  VYcenter);
fits.writeKey(fitmap, 'isPMRd', logical(isPMRd));
% Power meter level
fits.writeKey(fitmap, 'isAutScl', logical(isAutoScale));
ind = strfind(GnCubeNm ,'\');
fits.writeKey(fitmap, 'DatCube ', GnCubeNm(ind(end-1)+1:end));
ind = strfind(SRCubeNm ,'\');
fits.writeKey(fitmap, 'SRCube ', SRCubeNm(ind(end-1)+1:end));
ind = strfind(NormCubeNm ,'\');
fits.writeKey(fitmap, 'NormCube ', NormCubeNm(ind(end-1)+1:end));

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,meas);
fits.closeFile(fitmap);

%-- SAVE GAIN SETTINGS
fitmap= fits.createFile(GnCubeNm);
fits.createImg(fitmap,'double',size(scales));
%header data
ind = strfind(datCubeNm,'\');
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'Sample');
fits.writeKey(fitmap, 'NAX2', '1=FMTO_Scale;2=Gain;3=Bias;4=PMVal;5=Time');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,scales);
fits.closeFile(fitmap);

%-- SAVE SEMI REDUCED DATA
fitmap= fits.createFile(SRCubeNm);
fits.createImg(fitmap,'double',size(measScl));
%header data
ind = strfind(datCubeNm,'\');
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'Sample');
fits.writeKey(fitmap, 'NAX2', 'NRead bias subtracted and gain corrected');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,measScl);
fits.closeFile(fitmap);

%-- SAVE SEMI REDUCED NORM'D DATA
fitmap= fits.createFile(NormCubeNm);
fits.createImg(fitmap,'double',size(measScl2));
%header data
ind = strfind(datCubeNm,'\');
fits.writeKey(fitmap, 'DatCube ', datCubeNm(ind(end-1)+1:end));
fits.writeKey(fitmap, 'NAX1', 'Sample');
fits.writeKey(fitmap, 'NAX2', 'NRead bias subtracted, gain corrected, and normalized by PM');

fits.setCompressionType(fitmap,'NOCOMPRESS');
fits.writeImg(fitmap,measScl2);
fits.closeFile(fitmap);

%% Print results:
fprintf('\n Min value in Femto:             %f', min(measScl));
fprintf('\n Min value in Femto - Normed :   %f', min(measScl2));
fprintf('\n Max value in Femto:             %f', max(measScl));
fprintf('\n Min value in PM:                %f', min(pmRead));
fprintf('\n Max value in PM:                %f\n', max(pmRead));
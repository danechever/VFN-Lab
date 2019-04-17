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
%}
% Programm to run the piezo's in the VFN experiment to see the donut

%% add the function necessary
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\ControlCode\'));

close all
clear all

%% General settings 

% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\101418_FTO1';

% Experiment name for figures
expNm = '10Don_FixedFemtoLast4';

%~~ ZABER STUFF 
% Flag to enable zaber motion (including centering)
isZab = true;
isVortScan = false;

% Vortex center of scan properties [in mm]
VXcenter = 12.271445;%12.301445;
VYcenter = 13.134083;%13.111583;
Vbacklash= 0.005;

% Vortex scan properties
%To include center values as a point in scan; use ODD number of points
% !!! THIS IS VY !!!
VYpoints = 3;% Number of Y points in scan
% !!! THIS IS VX !!!
VXpoints = VYpoints; % Number of X points in scan

% Vortex step params in [mm]
vStepRange = 0.015;   %Vortex will be scanned +/- this value
                    % ie. scan will be: 
                    % [VXcen-vStepRa to VXcen+vStepRa] w/ Vxpoi steps
%~~ END ZABER STUFF 

%~~ PIEZO STUFF 
% Fiber center of scan in Volts
Xcenter = 77.42;%82.92;
Ycenter = 91.90;%86.97;
Zcenter = 55+5;

% Fiber scan properties
%To include center values as a point in scan; use EVEN number of points
Xpoints = 30;% number of X points in scan (actual scan will have +1)
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle
Zpoints = 1; %number of focci taken

% Fiber step sizes in Volts
refStep   = 8.0; % refStep is the step size for a full PSF with a 10*10 grid
StepSize  = refStep/(Xpoints/10);
ZStepSize = 2.5;
backlash  = 10; % Kept backlash in to see if it did anything
%~~ END PIEZO STUFF 

%~~ POWER METER STUFF 
Nread=10; %power samples at a given locaiton
% Add delay for settling time of detector
% NOTE::: Chose arbitrary delay values for now
delLookUp   = [[5 1]; [1 0.5]; [0 0.3]];  % Look-up table for delay values
delInd      = find(delLookUp<=StepSize,1);  % index of delay to use
Delay       =delLookUp(delInd,2); 
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)
%~~ END POWER METER STUFF 

%% Zaber Setup
if isVortScan && ~isZab
    error('isVortScan is true but vortex motion is disabled with isZab')
end

if isZab
    VFN_setUpZabers; % Instantiate the zabers

    if exist('ax14', 'var')
        vortX = ax14;
        clear ax14
    end
    if exist('ax63', 'var')
        vortY = ax63;
        clear ax63
    end
end

%% MDT (piezo) Setup
%Create serial port and define communication properties
    %COM# changes between computers; must check in Device Manager
MDT     = serial('COM4', 'BaudRate', 115200, 'InputBufferSize', 1500, ...
    'Terminator', {'CR' 'LF/CR'});
fopen(MDT);           %Connect port to device

%% Power meter setup
% this code is copied from Yinzi's code

%Setup Femto power meter

s = daq.createSession('ni');
addAnalogInputChannel(s,'Dev1', 0, 'Voltage');

s.Rate = 10; % rate is the rate of acquisition per second
s.DurationInSeconds = 0.1; 

%starting conversion scale (in V/W)
pm_scale = 10^6; % put here the current scale used to read the power


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

%Preallocate data matrix
meas=zeros(length(distX),length(distY),Nread,length(distZ),length(distVX), length(distVY)); 
% measPos = zeros(Xpoints,Ypoints,Zpoints,3);

if isZab
	%remove zaber X backlash
    VFN_Zab_move(vortX, VXcenter-Vbacklash);
end
for a = 1:length(distVX)
    if isZab
        %move vortex to new X position
        fprintf('\nVortX pos: %f',VFN_Zab_move(vortX, distVX(a)));
        %remove zaber y backlash
        VFN_Zab_move(vortY, VYcenter-Vbacklash);
    end
    for b = 1:length(distVY)
        if isZab
            %move vortex to new Y position
            fprintf('  VortY pos: %f\n',VFN_Zab_move(vortY, distVY(b)));
        end
        % Remove backlash in Z axis
        DE2_MDTVol(MDT, distZ(1)-backlash, 'z', 0);
        for k=1:length(distZ)
            
            DE2_MDTVol(MDT, distZ(k), 'z', 0);
            
            DE2_MDTVol(MDT, distX(1)-backlash, 'x', 0);
            for j=1:Xpoints+1
                fprintf('Vort %i of %i | Foc %i of %i | Column %i of %i\n', ...
                        (a-1)*length(distVY)+b,length(distVY)*length(distVX),...
                        k, length(distZ), j, length(distX))
                
                DE2_MDTVol(MDT, distX(j), 'x', 0);
                
                DE2_MDTVol(MDT, distY(1)-backlash, 'y', 0);
                for i=1:Ypoints+1
                    
                    DE2_MDTVol(MDT, distY(i), 'y', 0);
                    read=zeros(Nread,1);
                    %readPos = zeros(Nread, 3); % nReadx3 matrix for holding zaber positions
                    %col: (x, y, z), rows (reads)
                    tic;
                    pause(Delay);
                    for n=1:Nread
                        
                        read(n)=s.inputSingleScan; % taking power
                        if read(n) >10
                            warning('Power exceeds Femto limit of 10')
                        end
                        %readPos(n, 1) = ax14.getposition;
                        %readPos(n, 2) = ax63.getposition;
                        %readPos(n, 3) = ax93.getposition;
                        
                    end
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
                    test=toc;
                    meas(j,i,:,k,a,b)=read;
                    %meas(j,i,k)=mean(read);
                    %measPos(j,i,k,:) = mreadPos(:);
                end
                
            end
        end
    end
end
%meas=meas/pm_scale;

%% close the connections to all devices
%Set all three axes to starting position
DE2_MDTVol(MDT, Xcenter-backlash, 'x', 0); %xC holds current xVoltage
DE2_MDTVol(MDT, Xcenter, 'x', 0); %xC holds current xVoltage
DE2_MDTVol(MDT, Ycenter-backlash, 'y', 0); %yC holds current yVoltage
DE2_MDTVol(MDT, Ycenter, 'y', 0); %yC holds current yVoltage
DE2_MDTVol(MDT, Zcenter-backlash, 'z', 0); %zC holds current zVoltage
DE2_MDTVol(MDT, Zcenter, 'z', 0); %zC holds current zVoltage
%DE2_MDTVol(MDT, Zstart-backlash, 'z', 0); %zC holds current zVoltage
%DE2_MDTVol(MDT, Zstart, 'z', 0); %zC holds current zVoltage

fclose(MDT);
delete(MDT);
clear MDT;

if isZab
    %-- Return vortex to center position
    VFN_Zab_move(vortX, VXcenter-Vbacklash);
    fprintf('\nVortX pos: %f',VFN_Zab_move(vortX, VXcenter));
    VFN_Zab_move(vortY, VYcenter-Vbacklash);
    fprintf('  Vort Y pos: %f\n',VFN_Zab_move(vortY, VYcenter));
    
    %-- Clean up
    VFN_cleanUpZabers;
    if exist('vortX', 'var')
        clear vortX
    end
    if exist('vortY', 'var')
        clear vortY
    end
end

%% Display of data
for a = 1:length(distVX)
    for b = 1:length(distVY)
        for k = 1:length(distZ)
            figure();
            colormap('gray');
            %Not sure why transpose and "axis xy" are needed for correct orient.
            imagesc((distX-Xcenter),(distY-Ycenter),transpose(mean(meas(:,:,:,k,a,b),3)))
            titStr = sprintf('PWR Thru Fib, Z= %0.2f V',distZ(k));
            if isZab
                % Append vortex position if vortex was moved
                titStr = [titStr sprintf(', vort=(%0.5f, %0.5f)',distVX(a),distVY(b))]; %#ok<*UNRCH>
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
% section uncommanted. This section will creat a plot of the max power
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
meas1 = mean(meas,3);
plot(meas1(:,round(Ypoints/2+1)))
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
fitmap= fits.createFile([svFld filesep expNm '_' tag '.fits']);
fits.createImg(fitmap,'double',[length(distY) length(distX) Nread length(distZ) length(distVX) length(distVY)]);
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


fits.setCompressionType(fitmap,'NOCOMPRESS');

fits.writeImg(fitmap,permute(meas,[2,1,3,4,5,6]));
fits.closeFile(fitmap);

%% Report results
% Take average of all the Nread
meas1 = mean(meas,3);
% Crop the matrix to middle (cropVal-2)/cropVal region to be in donut null for min
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
         '\n  Foc ind = %3i, Foc = %6.2f     V', ...
         '\n  VX ind  = %3i, VX  =  %0.6f mm', ...
         '\n  VY ind  = %3i, VY  =  %0.6f mm'], ...
         I1+xmin-1, distX(I1+xmin-1), I2+ymin-1, distY(I2+ymin-1), ...
         I4, distZ(I4), I5, distVX(I5), I6, distVY(I6));


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
fprintf('\nMax/Min in min frame:    %f', mxVal/mnVal)
fprintf('\nMax/Min ratio-dark:      %f\n', (mxVal-0.005)/(mnVal-0.005))

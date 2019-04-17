%{ 
 Piezo Scan
 Code for scanning the fiber in the focal plane using the Piezo actuators
 It implements some of the things learned in Nicolas's experiments

 Removed all movement limit stuff since the Pzt's can only move 14um so
    they are extemely unlikely to collide with anything

*** This version saves a [X, Y, Nread, Focus] dimensional cube in the fits
      file instead of saving just the mean of the Nreads. 

*** This version does not use the zabers to move the vortex. This is from
      before I modified the system to automate the zabers
      ** Thus, there is a newer version of the code which allows you to
        control the vortex position using the zabers.
%}
% Programm to run the piezo's in the VFN experiment to see the donut

%% add the function necessary
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\ControlCode\'));

close all
clear all

%% General settings 

% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\092118_PZT7';

% Experiment name for figures
expNm = '28Don_ZabReplace11';

% Center of scan in Volts
Xcenter = 82.3;
Ycenter = 77.6;
Zcenter = 55;
%Zstart  = 75;

%To have center of psf on center of a read, use even number of points
Xpoints = 20;% number of points in X direction
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle
Zpoints = 1; %number of focus taken

%Step sizes in Volts
refStep   = 5; % refStep is the step size for a full PSF with a 10*10 grid
StepSize  = refStep/(Xpoints/10);
ZStepSize = 2.5;
backlash  = 10; % Kept backlash in to see if it did anything


Nread=10; %power samples at a given locaiton
% Add delay for settling time of detector
% NOTE::: Chose arbitrary delay values for now
delLookUp   = [[5 1]; [1 0.5]; [0 0.3]];  % Look-up table for delay values
delInd      = find(delLookUp<=StepSize,1);  % index of delay to use
Delay       =delLookUp(delInd,2); 
Delay = 0.1;
fprintf('Delay in use: %0.1f\n', Delay)

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
meas=zeros(Xpoints+1,Ypoints+1, Nread ,length(distZ)); % list where results will be stocked
% measPos = zeros(Xpoints,Ypoints,Zpoints,3);

% Remove backlash in Z axis
DE2_MDTVol(MDT, distZ(1)-backlash, 'z', 0); 

for k=1:length(distZ)
    
    DE2_MDTVol(MDT, distZ(k), 'z', 0);
    
    DE2_MDTVol(MDT, distX(1)-backlash, 'x', 0);
    for j=1:Xpoints+1
        fprintf('Foc %i of %i | Column %i of %i\n', k, length(distZ), j, Xpoints+1)
        
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
            meas(j,i,:,k)=read;
            %meas(j,i,k)=mean(read);
            %measPos(j,i,k,:) = mreadPos(:);
        end
        
    end
end
%meas=meas/pm_scale;

%% close the port
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

%% Display of data
for k = 1:length(distZ)
    figure(k);
    colormap('gray');
    %Not sure why transpose and "axis xy" are needed for correct orient.
    imagesc((distX-Xcenter),(distY-Ycenter),transpose(mean(meas(:,:,:,k),3)))
    title(sprintf('Power Through Fiber, Z = %0.5f V',distZ(k)))
    axis xy;axis square;
    cbr = colorbar;
    ylabel(cbr, 'Power [V]')
    set(gca,'TickDir','out')
    xlabel('Fiber X position [V]')
    ylabel('Fiber Y position [V]')
    tag = sprintf('%ix%i',Xpoints,Ypoints);
    if k ~= 1
        tag = sprintf([tag '_Foc%i'],k);
    end
    saveas(gcf, [svFld '\' expNm '_' tag '_CoupMap.png'])
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
fitmap= fits.createFile([svFld '\' expNm '_' tag '.fits']);
fits.createImg(fitmap,'double',[Xpoints+1 Ypoints+1 Nread length(distZ)]);
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

fits.setCompressionType(fitmap,'NOCOMPRESS');

fits.writeImg(fitmap,permute(meas,[2,1,3,4]));
fits.closeFile(fitmap);

%% Mins and Maxs
meas1 = mean(meas,3);
mnVal= min(min(min(meas1(Xpoints/5:4*Xpoints/5+1,Ypoints/5:4*Ypoints/5+1,:))));
fprintf('\nMin in center region: %f',mnVal)
mxVal = max(max(max(meas1)));
fprintf('\nMax:                  %f',mxVal)
ind = floor(find(meas1>= mxVal)/(Xpoints*Ypoints));
fprintf('\nFrame for max:        %d',ind)
fprintf('\nMax to Min ratio:     %f\n', mxVal/mnVal)
fprintf('\nMax/Min ratio-dark:   %f\n', (mxVal-0.005)/(mnVal-0.005))

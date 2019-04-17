%{ 
 Piezo Scan
 Code for scanning the fiber in the focal plane using the Piezo actuators
 It implements some of the things learned in Nicolas's experiments

 Removed all movement limit stuff since the Pzt's can only move 14um so
    they are extemely unlikely to collide with anything
%}
% Programm to run the piezo's in the VFN experiment to see the donut

%% add the function necessary
addpath(genpath('C:\Users\AOlab1\Desktop\DE2\VFN\PiezoCode')); % functions to control the zabers

close all
clear all

%% General settings 
%

% seems to move of a few micrometers between 2 days, if you loose it, give
% some defocus first then scan X and Y
% functions for controlling the zabers are in mm, these positions will also be in mm

% Directory for saving data:
svFld = 'C:\Users\AOlab1\Desktop\DE2\VFN\PupilVFNCoupling\090418_PZT';

% Experiment name for figures
expNm = '9PztScan_PSF';

% Center of scan in Volts
Xcenter = 75;
Ycenter = 97.7;
Zstart  = 83;

%To have center of psf on center of a read, use even number of points
Xpoints = 20;% number of points in X direction
Ypoints = Xpoints; % Ycenter/Xcenter will be an extra point in the middle
Zpoints = 1; %number of focus taken

%Step sizes in Volts
refStep   = 5; % refStep is the step size for a full PSF with a 10*10 grid
StepSize  = refStep/(Xpoints/10);
ZStepSize = 5;
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
pm_scale = 10^5; % put here the current scale used to read the power


%% Perform scan

meas=zeros(Xpoints,Ypoints,Zpoints); % list where results will be stocked
% measPos = zeros(Xpoints,Ypoints,Zpoints,3);

distX=(Xcenter-StepSize*Xpoints/2):StepSize:(Xcenter+StepSize*Xpoints/2);
distY=(Ycenter-StepSize*Ypoints/2):StepSize:(Ycenter+StepSize*Ypoints/2); 
distZ=Zstart:ZStepSize:(Zstart+ZStepSize*Zpoints);

% Remove backlash in Z axis
DE2_MDTVol(MDT, distZ(1)-backlash, 'z', 0); 

for k=1:Zpoints
    
    DE2_MDTVol(MDT, distZ(k), 'z', 0);
    
    DE2_MDTVol(MDT, distY(1)-backlash, 'y', 0);
    for j=1:Ypoints+1
        fprintf('Iteration #%i\n', j)
        
        DE2_MDTVol(MDT, distY(j), 'y', 0);
        
        DE2_MDTVol(MDT, distX(1)-backlash, 'x', 0);
        for i=1:Xpoints+1

            DE2_MDTVol(MDT, distX(i), 'x', 0);
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
            meas(i,j,k)=mean(read);
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
DE2_MDTVol(MDT, Zstart-backlash, 'z', 0); %zC holds current zVoltage
DE2_MDTVol(MDT, Zstart, 'z', 0); %zC holds current zVoltage

fclose(MDT);
delete(MDT);
clear MDT;

%% Display of data
for k = 1:Zpoints
    figure(k);
    colormap('gray');
    imagesc((distX-Xcenter),(distY-Ycenter),transpose((meas(:,:,k)))) % I don't know why I need the transpose here to have the good image, but it works that way, otherwise X and Y are swapped
    title(sprintf('Power Through Fiber, Z = %0.5f V',distZ(k)))
    axis xy;axis square;
    cbr = colorbar;
    ylabel(cbr, 'Power [V]')
    set(gca,'TickDir','out')
    xlabel('Fiber X position [V]')
    ylabel('Fiber Y position [V]')
    saveas(gcf, [svFld '\' expNm '_' sprintf('%ix%i_CoupMap.png',Xpoints,Ypoints)])
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
plot(meas(:,Ypoints+1))
xlabel('horizontal position of the fiber in the  focal plane in micrommeter')
ylabel('power in mWatts')
title(strcat('scan horizontal at ',num2str(Delay),'sec of delay'))
%% fits save
import matlab.io.*
fitmap= fits.createFile([svFld '\' expNm '_' sprintf('%ix%i.fits',Xpoints,Ypoints)]);
fits.createImg(fitmap,'double',[Xpoints+1 Ypoints+1 Zpoints]);
%header data
fits.writeKey(fitmap,'Xcenter ',Xcenter);
fits.writeKey(fitmap, 'Xpoints ',Xpoints);
fits.writeKey(fitmap,'Ycenter ',Ycenter);
fits.writeKey(fitmap,' Ypoints ',Ypoints);
fits.writeKey(fitmap,'Zstart ',Zstart);
fits.writeKey(fitmap,' Zpoints ',Zpoints);
fits.writeKey(fitmap,'XYsteps ',StepSize);
fits.writeKey(fitmap,' Zsteps ',ZStepSize);

fits.setCompressionType(fitmap,'NOCOMPRESS');
for k=1:Zpoints
fits.writeImg(fitmap,transpose(meas(:,:,k)));
end
fits.closeFile(fitmap);
%{
Auto Optimizer for Finding the Deepest Donut Null: 

******************************************************
- Arguments:
    NONE        USER MUST SET VARIOUS PARAMETERS:
      fldir     Folder to save data in
      lgnm      Update date at end of string
      h0:q      Initial conditions on DM shape
      UB:LB     Upper and Lower bounds on optimization region
- Returns:
    NONE        Outputs a .FITS cube with properties:
                    WidthxHeightxImage
                    200  x 200  x (heights -> phase)
- Dependencies:
    OPEN_multiDM()      Create DM connection
    UPDATE_multiDM()    Place new DM Map
    CLOSE_multiDM()     Disconnect DM
    uc480DotNet.dll     Library of CCD Camera Communications
    DE_DMW_Sin          Create Sinusoidal Map for DM
    DE_SupMinimizer_2   Iteration code around which to optimize
******************************************************


Author:    Daniel Echeverri
Last Modified:  09/07/2018
%}

STILL_WORKING_ON_MAKING_THIS_CODE

close all
clear all

global icount Pdat Idat ParT sb newp cam drv_inf img Xcr Ycr Samp

%% DEFINE OUTPUT DIRECTORY AND FILE; DEFINE SETUP PARAMETERS
%Output Directory Name: CHANGE LAST SECTION
fldDir  = 'C:\Users\COO_user\Documents\MATLAB\SpeckleNullingTests\2LDSpeckleResults0826_1f4';
%Set directory to folder with cubes for analysis
cd(fldDir);

%Name of TXT file for PM values
PMVnm   = 'PM_Values.txt';

%Name of TXT file for Raw Results
Rawnm   = 'RawResults.txt';

%Name of FITS file for Optimum Image
OptImnm = 'OptImage.fits';

%Name of Fits file for Image cube
Idatnm  = 'ImageTrials.fits';

%Number of Samples for PM data set; MUST BE AN INTEGER
Samp    = 500;

%Exposure time for camera
CExp    = .5;

%% Zaber Setup
if isVortScan && ~isZab
    error('isVortScan is true but vortex motion is disabled with isZab')
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

%if isAutoScale
%    % Take a sample reading at current power
%    readVal = mean(startForeground(s));
%    % Modify gain accordingly
%    [FMTO_scale, s] = VFN_FMTO_setAutoGain(s, readVal, FMTO_scale);
%else
    s = VFN_FMTO_setGain(s, FMTO_scale);
%end

fprintf('Current Gain setting: %i\n', FMTO_scale)

%% Red Thorlabs PM setup

%-- Connect to red PM only if it will be used
if isPMNorm
    % Find a VISA-USB object.
        % NOTE: check the SN on the device. The 2 student-lab PM100D's are:
        %   P0015519    and     P0015560
    obj1 = instrfind('Type', 'visa-usb', 'RsrcName', 'USB0::0x1313::0x8078::P0015519::0::INSTR', 'Tag', '');

    % Create the VISA-USB object if it does not exist
    % otherwise use the object that was found.
    if isempty(obj1)
        % NOTE: check the SN on the device. The 2 student-lab PM100D's are:
            %   P0015519    and     P0015560
        obj1 = visa('NI', 'USB0::0x1313::0x8078::P0015519::0::INSTR');
    else
        fclose(obj1);
        obj1 = obj1(1);
    end

    % Connect to instrument object, obj1.
    fopen(obj1);
end

%% Define Initial Parameters
Nactacross = 10;    %Number actuators across any given length
                        %VALUE IS AN ESTIMATE; NEED TO CONFIRM

% Set initial params: 
fibX = 75;          % Fiber X position
fibY = 75;          % Fiber Y position
fibZ = 75;          % Fiber Z position
vorX = 12.215000;   % Vortex X posiiton
vorY = 13.075000;   % Vortex Y position

% Vectorize initial params
par0 = [fibX, fibY, fibZ, vorX, vorY];

% SET NUMBER OF ITERATIONS TO PERFORM HERE:
itr = 3e3;

% Set optimization properties
fibBnd = 25;                                %Max distance from fibX/Y
focBnd = 50;                                %Max distance from fibZ
vorBnd = 0.25;                              %Max distance from focus
LB = [fibX-fibBnd,                          %Lower bounds
      fibY-fibBnd, 
      fibZ-focBnd, 
      vorX-vorBnd,
      vorY-vorBnd];
UB = [fibX+fibBnd,                          %Upper bounds
      fibY+fibBnd,
      fibZ+focBnd,
      vorX+vorBnd,
      vorY+vorBnd];
sTol = .0000001;                     %Smallest step size
fTol = 3e-9;                   %Stops iterating when delPower < fTol

%% Take Flat Data
%Create StringBuilder to retrieve and hold data; Set cap. to expected size
sb = System.Text.StringBuilder();
sb.Capacity = 4100;

%Create Matrix to hold images: Width(200)xHeight(200)x#Images            
    %Images have been cropped to reduce cube size
Idat = zeros(200,200,itr);    %Raw cropped images
%Create column vector to hold PM data
Pdat = zeros(itr, 1);
%Create Matrix to hold trial parameters
ParT = zeros(itr,4);

%Subtract 20 from Pdat for easy identification of empty elements later
    %In case fminsearch converges <itr iterations. Do not  want to process
    %or save excess images later
Pdat        = Pdat - 20;

%Reflatten cube for Flat image comparison
DE_DMW_Flat(0,drv_inf);
%Pause to ensure DM is flat; Note: DM mech response time <100us
pause(.05);

%Clear PM data store of any old data
newp.Write(4, 'PM:DS:CLear');

%_______Take Flat Data___________
%Capture data of flattened DM
%Acquire image
cam.Acquisition.Freeze(true);
%Take Data; Enable state returns to 0 automaticaly when 50 sample done
newp.Write(4, 'PM:DS:ENable 1');
pause(.05+Samp/10000);        %Collection takes Samp*.0001 seconds

%_______Process and Save Image Data______
%Extract image from RAM to array in Matlab
[ErrChk, tmp] = cam.Memory.CopyToArray(img.ID);
%Reshape array into plottable matrix
img.Data = reshape(uint8(tmp), [img.Width, img.Height, img.Bits/8]);
%Crop data to (200x200) around user input center; Save image
Idat(:,:,1)  = img.Data(int16(Ycr-99):int16(Ycr+100),int16(Xcr-99):int16(Xcr+100));
Iflat        = Idat(:,:,1);
%Draw cropped flat image
set(himg, 'CData', Idat(:,:,1));
axis tight equal off;
drawnow;

%_______Process and Save PM Data_________
sb.Clear();         %Ensure buffer is empty
%Place data in PM internal readable storage
newp.Write(4, ['PM:DS:Get? +' num2str(Samp)]);
pause(.05);         %Wait for data to finish storing in PM
%Loop to read arbitrary length data; Reads until no more data in PM
fl = 1;             %Flag to end reading           
str = '';           %variable to store data
while fl            %Loops until flag is tripped
    newp.ReadBinary(4, sb);             %Read 4096 characters from PM
    str = [str char(sb.ToString())];    %Append new read values
    %Change flag when the 'end of data' marker is found
    fl = ~strcmpi(str(end-12:end), ['end of data' 13 10]);
end
%Remove Header and footer from Data
str = str(281:end-13);        
%Save first PM value as flat; Take mean of 5000 points since signal
        %is noisy without the Analog and Digital Filter
Pdat(1) = mean(str2num(str));
Pflat   = Pdat(1);

%% Iterate, Optimize, and Process Data
%Define a counter to allow matching variables to each other
icount = 1;     %Start counter at 1 to compensate for flat image

%Use fminsearch to iterate and optimize
[Par, Pmin] = patternsearch(@DE_supMinimizer_2, par0, [], [], [], [], ... 
    LB, UB);%, [], optimoptions('fmincon', 'Algorithm', 'active-set'));%'MaxIterations', itr, 'FunctionTolerance', fTol));%,'StepTolerance', sTol));

%Remove empty parts of Image and Power arrays
tp   = find(Pdat == -20, 1) - 1;
Idat = Idat(:,:,1:tp);
ParT = ParT(1:tp, :);
Pdat = Pdat(1:tp);

%Get index of minimized value
ind = find(Pdat == Pmin, 1);

%Calculate Resulting Power Specs
OptP    = Pmin;
OptPSup = Pflat/Pmin;

%Match image to optimum power, show, and save as fits
OptIm   = Idat(:,:,ind);
figure()
imshow(OptIm);
caxis([0 255]);
fitswrite(OptIm, OptImnm);

%Save all images as single cube
fitswrite(Idat, Idatnm);

%Save PM values as vector
PMtxt = fopen(PMVnm, 'wt');
fprintf(PMtxt, '%d\n', Pdat);
fclose(PMtxt);

%_______Save Raw Results___________
REtxt = fopen(Rawnm, 'wt');
fprintf(REtxt, 'Folder with Results: %s \n', fldDir);
fprintf(REtxt, 'Number of Iterations: %d \n \n \n', tp-1);
fprintf(REtxt, 'Flat Power: %05.2e \n', Pflat);
%Save results for each iteration
for i = 2: icount
    fprintf(REtxt, 'Iteration: %d     ', i);
    fprintf(REtxt, 'Parameters: Amp. = %09.5f, Angle = %09.5f, Act. = %09.5f, Phase = %09.5f \n', ...
        ParT(i,1), ParT(i, 2), ParT(i,3), ParT(i,4));
    fprintf(REtxt, '       Power: %05.2e      Suppression Ratio: %05.2f \n \n', Pdat(i), Pflat/Pdat(i));
end
fprintf(REtxt, '\n \nOptimum Values: Iteration#%d       ', ind);
fprintf(REtxt, 'Parameters: Amp. = %09.5f, Angle = %09.5f, Act. = %09.5f, Phase = %09.5f \n', ...
        ParT(ind,1), ParT(ind, 2), ParT(ind,3), ParT(ind,4));
fprintf(REtxt, 'Power: %05.2e       Suppression Ratio: %05.2f', OptP, OptPSup);
fclose(REtxt);

%% Disconnect all Devices
%Return to MATLAB Directory
cd('C:\Users\COO_user\Documents\MATLAB');

%Disable and close Multi-DM driver USB connection
error_code = CLOSE_multiDM(drv_inf);
if error_code == 0
    fprintf('Successfully Closed DM\n');
else
    fprintf('DM CLOSE FAILED\n');
end

%Close PM
newp.CloseDevices();

%Close camera
if ~strcmp(char(cam.Exit), 'SUCCESS')
    fprintf('Could not close camera');
else
    fprintf('Successfully Close CCD\n');
end
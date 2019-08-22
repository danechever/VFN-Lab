% Varia_Pow_Read - Reads the power off of the red pm at the wavelengths and
% bandwidths indicated, has option for automation of power meter

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

% When the scan is done, 2 formatted matrices will be output in the command
% window, the Normalization values with the red pm to femto conversion and
% one without - the matrix without the conversion can be copied directly into the
% normalization matrix of the WavelengthScan.m to conduct a scan


%% Scan Parameters
% Wavelength centers being scanned
WLPoints = 11;
% Minimum WL and Maximum WL being scanned [Min,Max]
WLMinMax = [620,820];
% Max bandwidth of each center
MaxBW = 3;
% Total bandwidths scanned per center
BWPoints = 1;

%-- General Varia Settings
varPWR = 80;
varEMS = true;

% isPMNorm initiates and uses actuator for power meter
isPMNorm = true;

%% Initialize Varia
% Connect and instantiate NKT
VFN_setUpNKT;

%-- Set varia to desired initial setting
% Set power
VFN_NKT_setPowerLevel(NKT, varPWR);
% Set emission
VFN_NKT_setEmission(NKT, varEMS);

%% Connect to redPM
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
%% Read Values
% Checking if lights in room are off
fprintf(obj1, 'sense:correction:wavelength 780');
pmRead = nan(pmNread,1);
for ii = 1:pmNread
    pmRead(ii)=str2num(query(obj1, 'measure:power?'));
end
if(mean(pmRead) < 5e-9 || isPMNorm == false)
    % Move pm into beam
    VFN_Zab_move(pmX, 0);
    %Number of values being averaged
    pmNread = 100;
    % Wavelength centers being scanned
    WLCenters = linspace(WLMinMax(1),WLMinMax(2), WLPoints);
    %Bandwidths being measured for each center
    BWs = linspace(MaxBW/BWPoints,MaxBW,BWPoints);
    readVals = nan(length(WLCenters),length(BWs));
    rawRead = nan(length(WLCenters),length(BWs));
    for a = 1:length(WLCenters)
        for b = 1:BWPoints
            %-- Set Varia
            % Change bandwidth
            VFN_NKT_setWvlRange(NKT, WLCenters(a)-BWs(b)/2,WLCenters(a)+BWs(b)/2);
            % Change wavelength on red pm
            fprintf(obj1, ['sense:correction:wavelength ' num2str(WLCenters(a))]);
            disp(['Red PM changed to ' query(obj1, 'sense:correction:wavelength?')]);
            %pause for 60 seconds to let Varia stabilize
            pause(60)
            fprintf('Center: %3f Bandwidth: %3d\n',WLCenters(a), BWs(b));


            %Measuring norm value at single BW
            pmRead = nan(pmNread,1);
            for ii = 1:pmNread
                pmRead(ii)=str2num(query(obj1, 'measure:power?'));
            end
            pmRead1 = mean(pmRead(:))*0.997;% Updated on 2/8/19 based on 10/26/18
            pmRead1 = pmRead1*10^6;
            rawRead(a,b) = pmRead1;
            % Translate the calibration PM value from uW to V (at 10^6)
            pmRead1 = pmRead1*VFN_getFmtoResponsivity(WLCenters(a));%*.9;%0.63;
            fprintf('pmRead1: %6.5f    Raw: %6f \n', pmRead1, mean(pmRead))
            readVals(a,b) = pmRead1;
        end
    end
    disp(readVals')
else
    error('Room lights are on');
end
%% Report results
disp('pmNorm:');
fprintf('[');
for j = 1:(WLPoints-1)
    for i = 1:BWPoints
        fprintf(' %2.5f', readVals(j,i));
    end
    fprintf(';');
end
for i = 1:BWPoints
    fprintf(' %2.5f', readVals(WLPoints,i));
end
fprintf(']\n');

disp('pmNorm (raw):');
fprintf('[');
for j = 1:(WLPoints-1)
    for i = 1:BWPoints
        fprintf(' %2.5f', rawRead(j,i));
    end
    fprintf(';');
end
for i = 1:BWPoints
    fprintf(' %2.5f', rawRead(WLPoints,i));
end
fprintf(']\n');
%% Clean up NKT
VFN_cleanUpNKT;
%% Close redPM connection
fclose(obj1);
%% Clean up Zabers
if isPMNorm
    %-- Return vortex to center position
    %VFN_Zab_move(vortX, VXcenter-Vbacklash);
    %fprintf('\nVortX pos: %f',VFN_Zab_move(vortX, VXcenter));
    %VFN_Zab_move(vortY, VYcenter-Vbacklash);
    %fprintf('  Vort Y pos: %f',VFN_Zab_move(vortY, VYcenter));
    
    if isPMNorm
        %-- Move Calibration PM out of the beam
        fprintf('  |  PM Zab pos: %f\n', VFN_Zab_move(pmX, 25));
    else
        %fprintf('\n');  % terminate line from vortex recentering
        %-- Move fiber back to Zcenter
        %VFN_Zab_move(fibZ, Zcenter-Vbacklash);
        %VFN_Zab_move(fibZ, Zcenter);
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

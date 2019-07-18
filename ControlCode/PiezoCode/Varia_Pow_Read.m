%% Scan Parameters
% Wavelength centers being scanned
WLPoints = 1;
% Minimum WL and Maximum WL being scanned
WLMinMax = [635,635];
% Max bandwidth of each center
MaxBW = 60;
% Total bandwidths scanned per center
ScanPoi = 10;

%-- General Varia Settings
varPWR = 80;
varEMS = true;

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

%% Read Values
%Number of values being averaged
pmNread = 100;
% Wavelength centers being scanned
WLCenters = linspace(WLMinMax(1),WLMinMax(2), WLPoints);
%Bandwidths being measured for each center
BWs = linspace(MaxBW/ScanPoi,MaxBW,ScanPoi);
readVals = nan(length(WLCenters),length(BWs));
for a = 1:length(WLCenters)
    for b = 1:ScanPoi
        %-- Set Varia
        % Change bandwidth
        VFN_NKT_setWvlRange(NKT, WLCenters(a)-BWs(b)/2,WLCenters(a)+BWs(b)/2);
        %pause for 60 seconds to let Varia stabilize
        pause(60)
        disp(['Center: %3f Bandwidth: %3d',WLCenters(a), num2str(BWs(b))]);
        
        
        %Measuring norm value at single BW
        pmRead = nan(pmNread,1);
        for ii = 1:pmNread
            pmRead(ii)=str2num(query(obj1, 'measure:power?'));
        end
        pmRead1 = mean(pmRead(:))*0.997;% Updated on 2/8/19 based on 10/26/18
        pmRead1 = pmRead1*10^6;
        % Translate the calibration PM value from uW to V (at 10^6)
        pmRead1 = pmRead1*0.63;
        fprintf('pmRead1: %6.5f    Raw: %6f \n', pmRead1, mean(pmRead))
        readVals(a,b) = pmRead1;
    end
end
disp(readVals')

fprintf('[');
for i = 1:length(readVals)
    fprintf(' %2.5f', readVals(i));
end
fprintf('[');
%% Clean up NKT
VFN_cleanUpNKT;
%% Close redPM connection
fclose(obj1);
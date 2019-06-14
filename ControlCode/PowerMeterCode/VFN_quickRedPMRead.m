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

%% Read Value
pmNread = 100;
pmRead = nan(pmNread,1);
for ii = 1:pmNread
    pmRead(ii)=str2num(query(obj1, 'measure:power?'));
end
pmRead1 = mean(pmRead(:))*0.997;    % Updated on 2/8/19 based on 10/26/18
pmRead1 = pmRead1*10^6;
% Translate the calibration PM value from uW to V (at 10^6)
pmRead1 = pmRead1*0.63;
pmRead1
mean(pmRead)

%% Close connection
%%fclose(obj1)
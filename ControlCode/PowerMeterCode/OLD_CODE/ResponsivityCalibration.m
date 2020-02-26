%% Connect to redPM

% -Parameters
CalWvl = 780;

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

query(obj1, ['sense:correction:wavelength ' num2str(CalWvl)]);

VFN_setUpFMTO;  

s = VFN_FMTO_setGain(s, 5);
    
%% Read Value

pmNread = 500;
pmRead = nan(pmNread,1);
for ii = 1:pmNread
    pmRead(ii)=str2num(query(obj1, 'measure:power?'));
end
%pmRead1 = mean(pmRead(:))*0.997;    % Updated on 2/8/19 based on 10/26/18

disp(['pmRead: ' num2str(mean(pmRead))])
 
% Take a sample reading at current power
read = startForeground(s);
disp(['Femto: ' num2str(mean(read))])

%% Close connection

fclose(obj1)
%-- Clean up
FMTO_scale = 6;
s = VFN_FMTO_setGain(s, FMTO_scale);


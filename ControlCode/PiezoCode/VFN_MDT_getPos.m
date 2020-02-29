function rVol = VFN_MDT_getPos(MDT, Cha)
% VFN_MDT_GETPOS Function for checking the position of an MDT Piezo Axis
%   
%   - Assumes MDT object has already been created with VFN_setUpMDT
%   
%   * Based off the old DE2_MDTVol code
%
%   EXAMPLE:______
%   rVol = VFN_MDT_GETPOS(MDT, Cha)
%       MDT:        Serial object containing connection to MDT
%       Cha:        Axis to query ('x', 'y', or 'z') 
%       rVol:       Current position (in Volts)
%
%   See also VFN_setUpMDT, VFN_cleanUpMDT

% Check that Cha is valid
if ~any(strcmp(Cha, {'x','y','z','X','Y','Z'}))
    error('Provided axis is not valid, must be x, y or z')
end

% Set axis identifier to lower-case as needed by device
Cha = lower(Cha);

%Read Current Voltage
wri = [Cha 'voltage?'];  %Create read command; correct channel
fprintf(MDT, wri);              %Send read command
rVol = fscanf(MDT, '%s');       %Read result
rVol = str2double(rVol(3:end-1));%Convert to number

end
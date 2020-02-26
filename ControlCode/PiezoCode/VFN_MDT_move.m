function rVol = VFN_MDT_move(MDT, Vol, Cha)
% VFN_MDT_MOVE Function for moving an MDT Piezo Axis
%   
%   - Blocks execution until move completes
%   - Assumes MDT object has already been created with VFN_setUpMDT
%   - Throws error if desired position is beyond bounds of MDT.
%   * Allowable voltages: [0 to 150]
%   - Does large moves in 10V intervals to prevent fast, large voltage changes
% 
%   * Based off the old DE2_MDTVol code
%
%   EXAMPLE:______
%   rVol = VFN_MDT_MOVE(MDT, Vol, Cha)
%       MDT:        Serial object containing connection to MDT
%       Vol:        position to move to (in Volts)
%       Cha:        Axis to move ('x', 'y', or 'z') 
%       rVol:       resulting position (in Volts)
%
%   See also VFN_setUpMDT, VFN_cleanUpMDT

% Set axis identifier to lower-case as needed by device
Cha = lower(Cha);

%Read Current Voltage
wri = [Cha 'voltage?'];  %Create read command; correct channel
fprintf(MDT, wri);              %Send read command
rVol = fscanf(MDT, '%s');       %Read result
rVol = str2double(rVol(3:end-1));%Convert to number

% Check if desired voltage is within allowable region [0,150]
if (Vol < 0 || Vol >150)
    error('Desired MDT voltage is beyond bounds [0 - 150]');
end

%Loop to perform smooth change to desired voltage; 10V steps
fla  = true;                %Exit flag
lct  = 0;                   %Counter to avoid infinite loops
while fla
    %Check distance from current voltage to desired voltage
    dV   = Vol - rVol;
    %If too far, do 10V step, else do dV step
    if abs(dV) > 10
        Volw = rVol + 10*dV/abs(dV);%dV/abs(dV) defines direction
    else
        Volw = Vol;%rVol + dV;
    end

    %Create write command; correct channel and voltage with formatting
    wri = [Cha 'voltage=' num2str(Volw, '%-.2f')];
    fprintf(MDT, wri);              %Send message on serial port
    wri = [Cha 'voltage?'];         %Create read command
    pause(.075);                      %Wait for voltage to settle
    fprintf(MDT, wri);              %Send read command
    rVol = fscanf(MDT, '%s');       %Read result
    rVol = str2double(rVol(4:end-1));%Convert to number

    %Check whether further iteration is needed
    fla = abs(Vol-rVol)>.5;         %.5 to avoid infinite loop
                              %BC Change rom .3 for less precision

    %Increment and check counter; quite to avoid infinite loop
    lct = lct + 1;
    if lct == 30 %BC Change from 20 to give more flexibility
        fla = false;
        error('Infinite loop when writing to MDT\n');
    end
end

end
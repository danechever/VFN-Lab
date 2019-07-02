%{
Piezocontroller Communications: Voltage Control
- Sets or reads voltage to given channel on Pizeocontroller (MDT)
- Read simply reads the value from the MDT
- Write performs a write and then a read to return the value
* Writes are performed in 10V steps to avoid large voltage changes
- Only a RW=0 performs a write to avoid accidental writes; all else results
    in a read
- Automatically pushes writes into desired bounds [0, 150]
*** ASSUMES serial port already created and opened w/ correct properties
        MDT = serial(...);
        fopen(MDT);

******************************************************
- Arguments:
    ser         serial object to write to (MDT)
    Vol         New value for voltage       
    Cha         Channel to read/write from/to:
                 x = xvoltage
                 y = yvoltage
                 z = zvoltage
    RW          Integer to choose read vs. write:
                 0   => Write
                 !=0 => Read
                
- Returns:
    rVol        Value of channel written to:
                 in Read: Value from device
                 in Write: Value returned after writing
                 
- Dependencies:
    None        Uses only default MATLAB functions
******************************************************


Authors:    Daniel Echeverri and Michael Randolph
Last Modified:  09/04/2018
%}

function rVol = DE2_MDTVol(MDT, Vol, Cha, RW)
    Cha = lower(Cha);

    %Read Current Voltage
    wri = [Cha 'voltage?'];  %Create read command; correct channel
    fprintf(MDT, wri);              %Send read command
    rVol = fscanf(MDT, '%s');       %Read result
    rVol = str2double(rVol(3:end-1));%Convert to number
    
    %Perform a write if RW is 0; overwrites rVal with value after writing
    if RW == 0
        %Push voltage towards allowable region [0,150]
        if Vol < 0
            Vol = 0;
        end
        if Vol > 150
            Vol = 150;
        end
        
        %Loop to perform smooth change to desired voltage; 10V steps
        fla  = 1;                   %Exit flag
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
            fla = abs(Vol-rVol)>.3;         %.05 to avoid infinite loop
            
            %Increment and check counter; quite to avoid infinite loop
            lct = lct + 1;
            if lct == 20
                fla = 0;
                error('Infinite loop when writing to MDT\n');
            end
        end
    end
end
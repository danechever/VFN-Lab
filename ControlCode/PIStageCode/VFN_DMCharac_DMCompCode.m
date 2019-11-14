%{
Script showing basic serial communication between two computers. In this
case, it is structured for characterizing the PI-Q-545 stage with the Zygo 
but it can easily be adapted for other applications.

::: This is built on the code originally written for characterizing the DM

::: This version is the DM Comp's shell version

Has 3 sections
1) Basic functions (explains serial communication)
2) DM Computer Code (shell for DM computer code)
3) Zygo Computer Code (shell for zygo computer code)

The DM and Zygo code implement an artificial handshaking to ensure that
    data capture is synchronized. This is done by sending keywords from one
    to the other and enables a simple timeout implementation while ensuring
    a reliable way to exit the main loop.
Also implements a master-slave (primary-secondary) system with the DM 
    Comp as primary, controlling when data capture begins and ends on the
    Zygo Comp (secondary)

%NOTE:
    - The computers must be connected using a null-modem cable:
        ie. serial cable with pins 2 and 3 crossed on each end
    - If the computer does not have a serial port (DM comp) you'll need to 
        use a serial-usb converter. There are several in the lab; I 
        recommend using the long one with the clear, patterned wire. The 
        others do not work as well.

    - Using a serial port is very similar to writing to a txt file: you
        create a serial object, open the port, print your data, read you 
        data, and then close the port when you're done.
 
******************************************************
- Arguments:
    NONE            *You need to change certain parts within the code as
                        prompted by comments in the code below.

- Returns:
    NONE            *No formal returned values. Creates zygo measurement
                        data on the Zygo computer.
 
- Dependencies:
    JR_Shape_Measurement    Function to take zygo measurements  
        (and dependencies related to this function)
******************************************************

Initial Version:    Daniel Echeverri, 02/10/2018
Last Edit:          Daniel Echeverri
Last Modified:      02/10/2018
%}

% IF USING A SERIAL-USB ADAPTER, ADAPTER MUST BE CONNECTED AND RECOGNIZED
    % BEFORE STARTING MATLAB

%% Basic Serial communication structure
    %Basic functions. Actual suggested code is in the next sections
    %Feel free to jump to next sections; this one just explains what I'm
        %doing in the other ones
%{    
%Create serial port object. 
s = serial('COM15', 'BaudRate', 921600, 'Parity', 'none', 'StopBits', 1, 'FlowControl', 'none');
    %COM        = serial port. Determined from the device manager.
    %BaudRate   = data transfer rate (bits/second). 
    %Parity     = Method of confirming successful tx. Not needed here.
    %StopBits   = Framing bit. Keep as 1.
    %FlowControl= Handshaking type. Not needed here. In fact, with the wrong type of cable, this'll cause problems.
    
%Open the port
fopen(s);

%OPEN THE SERIAL PORT ON YOUR OTHER COMPUTER BEFORE CONTINUING

%Send (ASCII) Data:_________
%* Messages themselves are sent as ascii strings so you'll need to convert
        %before transmission
%1) Single string 
fprintf(s, 'Hello World from PC1');
%2) Number (use fprintf formatspecs - %f = float, %i = int, etc.)
fprintf(s, '%f', 3.7);
%3) Vector 
vect = [1; 2.3; 3.7; 4.0; -5];
fprintf(s, '%f', vect);
    %This will send the whole vector separated by the default terminator (see below)

%Read (ASCII) Data:________
%1) Read text
rd = fscanf(s);
%2) Read numbers
rdInt = fscanf(s, '%i');
%3) Read 3 lines
rdFlt = fscanf(s, '%f', 3);

%Display all settings
set(s)
%Display with current values
get(s)
%Change the timeout (amount of time system waits while reading a message)
s.Timeout = 3;  % (seconds)
%Change the terminator (character appended to the end of every message)
s.Terminator = CR\LF;    %carriage return\linefeed 
%Change the input buffer size (number of bytes that can be held stored by
    %the system at any given moment. fscanf actually reads and removes bytes
    %from this buffer)
s.InputBufferSize = 1024; %default is 512
%Check if there is data in the read buffer
nBytes = s.BytesAvailble;

%Close the port
fclose(s);
%}
%% Code for DM comp (primary)
%CHECK THAT WE CAN EXTEND THE ZYGO WAIT TIME IN METROPRO

%::::::::::::Prep Serial Comms::::::::::::
%Check if serial object already exists
if exist('s','var')
    fclose(s);
    delete(s);
    clear s;
end

%Start comms ***Max Baud for Zygo comp is 115200; both should be at same baud
s = serial('COM15', 'BaudRate', 115200, 'Parity', 'none', 'StopBits', 1, 'FlowControl', 'none');
fopen(s);
tmt     = 3000;     %Artificial timeout for my handshaking. 
                    %3000 = 300s = 5 min
%::::::::::::Serial Comms Prepped::::::::::::

                    
%::::::::::::Wait for Zygo Comp Connection::::::::::::
%My implementation of artificail handshaking and timeout
fprintf('\nWaiting for Zygo Computer Connection')
conFlg  = 0;        %flag to mark communication success
for i = 0:tmt       %iterate until timeout or message is received
    %send 'Start' every 5s in case Zygo comp serial port was not open b4
    if mod(i,50) == 0
        fprintf(s, 'Start');
        fprintf('\n   Waiting... itr %i', i)
    end
    if s.BytesAvailable ~= 0    %Check for received messages
        rd = fscanf(s, '%s');
        if strcmpi(rd,'Start')        %Exit if 'Start' is received
            conFlg = 1;
            fprintf(s, 'Start');    %Send final 'Start' in case other comp started in between 'Start' sends
            break
        end
    end
    pause(.1)
end
%Throw error to quit code if connection fails
if ~conFlg
    error('Zygo comp did not respond')
end
pause(.5)
%Clear input buffer of all data to ensure synchronization
flushinput(s)
pause(1)
%::::::::::::Zygo Comp Connection Done::::::::::::::


%*AT THIS POINT, BOTH COMPUTERS HAVE CONFIRMED COMMUNICATION*
fprintf('\nConnected to Zygo Computer')


%::::::::::::PIStage Prep Code::::::::::::
%RUN PIStage PREP CODE HERE (<tmt on Zygo Comp or else Zygo Comp will timeout)
VFN_PIStage_init;
VFN_PIStage_home;
%::::::::::::PIStage Prep Done::::::::::::


%::::::::::::Main Loop Start::::::::::::
%Create a variable with the number of iterations to perform (number of measurements)
msrs = 21;
mvs  = logspace(-5.3,-3,msrs-1);  %1e-3 (1micron) = a ~0.018Lambda shift
                                % Thus, (-5.3,-3) will cover 5nm-1micron
mvs  = [0 mvs];     % Add a 0-shift value to the list of positions
startPos = 0.3383; % Reference position (ie. 0-point for shifts)
for i = 1:msrs
    
    %::::::::::::PIStage code goes here::::::::::::
    PIpos=startPos-mvs(i);
    VFN_PIStage_move;
        
    %Create string for filename on Zygo comp (ex: using x,y coord of act and voltage)
        %This can be whatever you want it to be...
    flnm    = sprintf('msr_%i_Shift_-%06.1f',i, mvs(i)*1e6);
    %:::::::::::::PIStage poke stable::::::::::::::
       
    
    %::::::::::::Cue Zygo measure:::::::::::::
    %Send message marking start of measure (artificial handshaking)
    fprintf(s, 'Measure');
    pause(.5)
    %Send filename
    fprintf(s, flnm);
    fprintf(['\nWaiting on Zygo to measure: ' flnm])
    %Wait for Zygo to reply that it is done
        %Zygo replies with flnm to ensure that measures are synchronized
    conFlg  = 0;        %flag to mark communication success
    fprintf('\n   Waiting for measure to complete... ')
    for j = 1:tmt
        if mod(j,50) == 0
            fprintf('%i ', j/10)    
        end
        if s.BytesAvailable ~= 0    %Check for received messages
            rd = fscanf(s, '%s');
            if strcmpi(rd, flnm)  %Exit if correct filename is received
                conFlg = 1;
                break
            end
        end
        pause(.1)
    end
    %Throw error to quit code if connection fails
    if ~conFlg
        % Return stage to starting position
        PIpos=startPos; VFN_PIStage_move;
        % Cleanup stage connection
        VFN_PIStage_cleanUp;
        error('Zygo comp did not respond or responded incorrectly')
    end
    fprintf(['\n' flnm ' measurement done\n'])
    %::::::::::::Zygo measure Done::::::::::::
end
%::::::::::::Main Loop Done::::::::::::


%::::::::::::Tell Zygo Comp to Stop::::::::::::
fprintf(s, 'Stop');    %Needed in primary-secondary structure
fprintf('\nStop Sent\n')
%::::::::::::Zygo Comp Told Stop::::::::::::::


%::::::::::::Close Serial Comms::::::::::::
fclose(s);
delete(s);
clear s;
%::::::::::::Serial Comms Closed::::::::::::


%::::::::::::PIStage Close Code::::::::::::
% Return stage to starting position
PIpos=startPos; VFN_PIStage_move;
% Cleanup stage connection
VFN_PIStage_cleanUp;
%::::::::::::PIStage Close Done::::::::::::
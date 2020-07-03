%Set up the Thorlabs MDT for use
%
%- Connects to the MDT controller which runs the Thorlabs piezo actuators
%- Creates an "MDT" variable containing pertinent object information
%
% See also VFN_MDT_move, VFN_cleanUpMDT
%
%  author: D. Echeverri (dechever@caltech.edu)
%  created: 02/25/2020

% Create serial port and define communication properties
    %COM# changes between computers; must check in Device Manager
MDT     = serial('COM5', 'BaudRate', 115200, 'InputBufferSize', 1500, ...
    'Terminator', {'CR' 'LF/CR'}, 'Timeout', 1);

% Connect to device
fopen(MDT);
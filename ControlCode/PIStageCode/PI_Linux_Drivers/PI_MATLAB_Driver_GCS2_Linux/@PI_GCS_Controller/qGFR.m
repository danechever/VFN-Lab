function data = qGFR (c, tabelIds, startPoint, numberOfPoints)
% Read the data stored in the data recorder.
%
%   SYNTAX
%   data = PIdevice.qDRR (tabelIds, startPoint, numberOfPoints)
%
%   DESCRIPTION
%   Get the data store#d in the data recorder of the controller.
%
%   INPUT
%       tabelIds              IDs of the table whose data you want to read.
%       (integer array)
%
%       startPoint            Start point in the data recorder table,
%       (double)              starts with index 1.
%
%       numberOfPoints        Number of points to read. Set to "-1" to
%       (double)              read all points available (works for most
%                             controllers, but not all).
%
%   OUTPUT
%       data                  Data stored by the data recorder. The first
%                             row -> data (:, 1) will contain the time
%                             base, the following row(s) will contain the
%                             data recorded.
%
%   EXAMPLE
%       Get all values stored in table 1, 2 and 3.
%       data = PIdevice.qDRR ([1 2 3], 1, -1)
%
%       Get the values stored in table 1.
%       You will get values with index number 11, 12, 13, ... 59, 60.
%       data = PIdevice.qDRR (1, 11, 50)
%

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

functionName = [ c.libalias, '_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found', functionName), end;

% Create variables for C interface
piTables     = libpointer('int32Ptr',tabelIds);
nTables      = length(tabelIds);
headerLength = 10000;
header       = blanks ( headerLength + 1 );
ppdData      = libpointer ( 'doublePtr' );

try
    % call C interface
    [bRet,~,ppdData,header] ...
        = calllib(c.libalias,functionName,c.ID,piTables,nTables,startPoint,numberOfPoints,ppdData,header,headerLength);

    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end

%% Parse GCS-Array Header
headerlines = regexp ( header,'\n','split' );

sampletime   = -1;
timeColumnIndex = -1;

for n = 1:length(headerlines)
    currentline = headerlines{n};
    if(~isempty(strfind(currentline,'DIM')))
        nTables = str2num(currentline(strfind(currentline,'=')+1:end));
    end
    if(~isempty(strfind(currentline,'NDATA')))
        numberOfPoints = str2num(currentline(strfind(currentline,'=')+1:end));
    end
end

%% Parse GCS-Array Data

% Wait until all Data is received
numberOfDatapoints = 0;
while ( numberOfDatapoints < ( nTables * numberOfPoints ) )
    pause ( 0.1 );
    numberOfDatapoints =  GetAsyncBufferIndex ( c );
end

setdatatype ( ppdData,'doublePtr', nTables, numberOfPoints );
data = ppdData.Value';


function [ returnValues ] = qFSF ( c, szAxis )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
dForceValue1Array = 0;
pForceValue1Array = libpointer ( 'doublePtr', dForceValue1Array );
dPositionOffsetArray = 0;
pPositionOffsetArray = libpointer ( 'doublePtr', dPositionOffsetArray );
dForceValue2Array = 0;
pForceValue2Array = libpointer ( 'doublePtr', dForceValue2Array );

try
    % call C interface
    [ bRet, ~, dForceValue1Array, dPositionOffsetArray, dForceValue2Array ] = calllib ( c.libalias, functionName, c.ID, szAxis, pForceValue1Array, pPositionOffsetArray, pForceValue2Array );
    
    returnValues = { dForceValue1Array, dPositionOffsetArray, dForceValue2Array };
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
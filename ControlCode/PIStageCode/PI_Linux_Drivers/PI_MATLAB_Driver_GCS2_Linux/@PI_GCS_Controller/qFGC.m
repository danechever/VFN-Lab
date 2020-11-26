function [ dScanAxisCenterValueArray, dStepAxisCenterValueArray ]  = qFGC ( c, szProcessIds )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

len = GetNrAxesInString ( c, szProcessIds );
if(len == 0)
    return;
end

% Create variables for C interface
dScanAxisCenterValueArray = zeros ( len, 1 );
dStepAxisCenterValueArray = zeros ( len, 1 );

pdScanAxisCenterValueArray = libpointer ( 'doublePtr', dScanAxisCenterValueArray );
pdStepAxisCenterValueArray = libpointer ( 'doublePtr', dStepAxisCenterValueArray );


try
    % call C interface
    [ bRet, ~, dScanAxisCenterValueArray, dStepAxisCenterValueArray ] ...
        = calllib ( c.libalias, functionName, c.ID, szProcessIds, pdScanAxisCenterValueArray, pdStepAxisCenterValueArray );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
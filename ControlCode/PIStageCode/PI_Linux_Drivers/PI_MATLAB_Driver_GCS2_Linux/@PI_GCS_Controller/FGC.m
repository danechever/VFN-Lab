function FGC ( c, szProcessIds, dScanAxisCenterValueArray, dStepAxisCenterValueArray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

% Create variables for C interface
pdScanAxisCenterValueArray = libpointer('doublePtr', dScanAxisCenterValueArray);
pdStepAxisCenterValueArray = libpointer('doublePtr', dStepAxisCenterValueArray);

try
    % call C interface
    bRet = calllib ( c.libalias, functionName, c.ID, szProcessIds, pdScanAxisCenterValueArray, pdStepAxisCenterValueArray );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
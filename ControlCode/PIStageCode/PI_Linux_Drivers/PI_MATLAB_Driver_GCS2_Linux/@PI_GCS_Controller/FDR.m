function FDR ( c, szScanRoutineName, szScanAxis, dScanAxisRange, szStepAxis, dStepAxisRange, szParameters )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

try
    % call C interface
    bRet = calllib ( c.libalias, functionName, c.ID, szScanRoutineName, szScanAxis, dScanAxisRange, szStepAxis, dStepAxisRange, szParameters );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
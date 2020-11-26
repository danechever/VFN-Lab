function FRP ( c, szScanRoutineNames, iRoutineActionsArray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

% Create variables for C interface
piRoutineActionsArray = libpointer ( 'int32Ptr', iRoutineActionsArray );

try
    % call C interface
    [ bRet ] = calllib ( c.libalias, functionName, c.ID, szScanRoutineNames, piRoutineActionsArray );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
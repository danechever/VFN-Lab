function iRoutineActionsArray = qFRP ( c, szScanRoutineNames )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
len = GetNrAxesInString ( c, szScanRoutineNames );
if(len == 0)
    return;
end

iRoutineActionsArray = zeros ( len, 1 );
piRoutineActionsArray = libpointer ( 'int32Ptr', iRoutineActionsArray );


try
    % call C interface
    [ bRet, ~, iRoutineActionsArray ] = calllib ( c.libalias, functionName, c.ID, szScanRoutineNames, piRoutineActionsArray );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
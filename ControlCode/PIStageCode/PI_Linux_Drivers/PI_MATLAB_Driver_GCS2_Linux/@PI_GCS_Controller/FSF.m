function FSF ( c, szAxis, dForceValue1, dPositionOffset, iUseForceValue2, dForceValue2 )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface


try
    % call C interface
    [ bRet ] = calllib ( c.libalias, functionName, c.ID, szAxis, dForceValue1, dPositionOffset, iUseForceValue2, dForceValue2 );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
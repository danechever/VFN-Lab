function WCL ( c, iWaveTableIdsArray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

piWaveTableIdsArray = libpointer( 'int32Ptr', iWaveTableIdsArray );
iNumberOfParameters = length ( iWaveTableIdsArray );

try
    % call C interface
    bRet = calllib ( c.libalias, functionName, c.ID, piWaveTableIdsArray, iNumberOfParameters );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
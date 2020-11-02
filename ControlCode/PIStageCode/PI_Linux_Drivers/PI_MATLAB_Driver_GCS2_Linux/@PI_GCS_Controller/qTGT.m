function iTrajectoryTiming = qTGT ( c )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
iTrajectoryTiming = 0;
piTrajectoryTiming = libpointer ( 'int32Ptr', iTrajectoryTiming );

try
    % call C interface
    [ bRet, iTrajectoryTiming ] = calllib ( c.libalias, functionName, c.ID, piTrajectoryTiming );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
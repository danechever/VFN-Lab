function NumberOfRecordedFrequencyMeasurements = qGFL ( c )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
NumberOfRecordedFrequencyMeasurements = 0;
piRoutineActionsArray = libpointer ( 'int32Ptr', NumberOfRecordedFrequencyMeasurements );

try
    % call C interface
    [ bRet, NumberOfRecordedFrequencyMeasurements ] = calllib ( c.libalias, functionName, c.ID, piRoutineActionsArray );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
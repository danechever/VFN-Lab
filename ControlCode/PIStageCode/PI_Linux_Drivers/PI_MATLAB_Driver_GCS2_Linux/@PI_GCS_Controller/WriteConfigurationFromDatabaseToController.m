function WriteConfigurationFromDatabaseToController ( c, szFilter, szConfigurationName )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
bufferSize = 20000;
szWarnings = blanks(20000);


% try
% call C interface
[ returnValue, ~, ~, warnings ] ...
    = calllib ( c.libalias, functionName, c.ID ...
    , szFilter, szConfigurationName, szWarnings, bufferSize );

if ( 0 == returnValue )
    gcsError.Code = GetError ( c );
    gcsError.Description = TranslateError ( c, gcsError.Code );
    
    completeString = sprintf('Error Code: %d\n', gcsError.Code);
    completeString = sprintf('%s%s\n\n', completeString, gcsError.Description);
    completeString = sprintf('%s%s', completeString, warnings);
    
    if (  gcsError.Code == -10013 ...   % PI_PARAMETER_DB_AND_HPA_MISMATCH_LOOSE
       || gcsError.Code == -10014 )     % PI_PARAMETER_DB_FAILED_TO_SET_PARAMETERS_CORRECTLY
        warning (completeString);        
    else
        error (completeString);
    end
end

% catch ME

% end
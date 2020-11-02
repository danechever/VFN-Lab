function szConfigurations = GetAvailableControllerConfigurationsFromDatabase ( c )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
bufferSize = 20000;
szConfigurations = blanks(20000);


try
    % call C interface
    [ bRet, szConfigurations ] = calllib ( c.libalias, functionName, c.ID , szConfigurations, bufferSize );
   
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error (ME.message);
end
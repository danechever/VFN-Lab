function [ szAnswer ] = qFRC ( c, szProcessIdBase )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

% Create variables for C interface
szAnswer = blanks ( 10001 );

try
    % call C interface
    [ bRet, ~, szAnswer ] = calllib ( c.libalias, functionName, c.ID, szProcessIdBase, szAnswer, 10000 );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
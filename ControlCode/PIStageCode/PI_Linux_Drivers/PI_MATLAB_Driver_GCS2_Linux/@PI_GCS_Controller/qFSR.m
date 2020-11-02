function [ bValueArray ] = qFSR ( c, szAxis )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
bValueArray = 0;
pbValueArray = libpointer ( 'int32Ptr', bValueArray );


try
    % call C interface
    [ bRet, ~, bValueArray ] = calllib ( c.libalias, functionName, c.ID, szAxis, pbValueArray );
   
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
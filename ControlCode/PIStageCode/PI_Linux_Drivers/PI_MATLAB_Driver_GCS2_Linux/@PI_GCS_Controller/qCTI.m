function szValueArray = qCTI ( c, iTriggerInputIds, iTriggerParameters )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piTriggerInputIds = libpointer('int32Ptr', iTriggerInputIds);
piTriggerParameters = libpointer('int32Ptr', iTriggerParameters);
iArraySize = length ( iTriggerInputIds );
szValueArray = blanks ( 10001 );


try
    % call C interface
    [ bRet, ~, ~, szValueArray ] = calllib ( c.libalias, functionName, c.ID, ...
        piTriggerInputIds, piTriggerParameters, szValueArray, iArraySize, 10000 );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
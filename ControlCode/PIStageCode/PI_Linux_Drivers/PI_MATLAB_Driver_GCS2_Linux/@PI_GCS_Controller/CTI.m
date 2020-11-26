function CTI ( c, iTriggerInputIds, iTriggerParameterArray, szValueArray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piTriggerInputIds = libpointer('int32Ptr', iTriggerInputIds);
piTriggerParameterArray = libpointer('int32Ptr', iTriggerParameterArray);
iArraySize = length(iTriggerParameterArray);


try
    % call C interface
    [ bRet ] = calllib ( c.libalias, functionName, c.ID, piTriggerInputIds, piTriggerParameterArray, szValueArray, iArraySize );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
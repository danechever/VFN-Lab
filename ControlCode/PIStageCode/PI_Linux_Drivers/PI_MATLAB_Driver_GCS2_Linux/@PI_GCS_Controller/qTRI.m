function bTriggerState = qTRI ( c, iTriggerInputIds )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piTriggerInputIds = libpointer('int32Ptr', iTriggerInputIds);
bTriggerState = zeros ( size ( iTriggerInputIds ) );
pbTriggerState = libpointer('int32Ptr', bTriggerState);
iNumberOfInputIds = length ( iTriggerInputIds );


try
    % call C interface
    [ bRet, ~, bTriggerState ] = calllib ( c.libalias, functionName, c.ID, piTriggerInputIds, pbTriggerState, iNumberOfInputIds );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
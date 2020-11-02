function dCalculatedInputValueArray = qTCI ( c, iFastAlignmentInputIdsArray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piFastAlignmentInputIdsArray  = libpointer ( 'int32Ptr', iFastAlignmentInputIdsArray );

iArraySize = length ( iFastAlignmentInputIdsArray );
dCalculatedInputValueArray = zeros ( iArraySize, 1 );
pdCalculatedInputValueArray = libpointer('doublePtr', dCalculatedInputValueArray);


try
    % call C interface
    [ bRet, ~, dCalculatedInputValueArray ] = calllib ( c.libalias, functionName, c.ID, piFastAlignmentInputIdsArray, pdCalculatedInputValueArray, iArraySize );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
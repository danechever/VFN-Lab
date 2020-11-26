function szAnswer = qSIC ( c, iFastAlignmentInputIdsArray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piFastAlignmentInputIdsArray = libpointer('int32Ptr', iFastAlignmentInputIdsArray);
iNumberOfInputIds = length ( iFastAlignmentInputIdsArray );
szAnswer = blanks ( 10001 );

try
    % call C interface
    [ bRet, ~, szAnswer ] = calllib ( c.libalias, functionName, c.ID, piFastAlignmentInputIdsArray, iNumberOfInputIds, szAnswer, 10000 );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
function SIC ( c, iFastAlignmentInputId, iCalcType, dParameters )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;

% Create variables for C interface
pdParameters = libpointer('doublePtr', dParameters);
iNumberOfParameters = length ( dParameters );

try
    % call C interface
    [ bRet ] = calllib ( c.libalias, functionName, c.ID, iFastAlignmentInputId, iCalcType, pdParameters, iNumberOfParameters );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
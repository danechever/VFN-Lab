function TGA ( c, iTrajectoriesArray, dValarray )

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piTrajectoriesArray = libpointer('int32Ptr', iTrajectoriesArray);
pdValarray = libpointer('doublePtr', dValarray);
iArraySize = length(dValarray);


try
    % call C interface
    [ bRet ] = calllib ( c.libalias, functionName, c.ID, piTrajectoriesArray, pdValarray, iArraySize );
    
    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end
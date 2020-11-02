function [bRet, charBuffer] = qKEN(c, szNamesOfCoordSystems)
% function [bRet, szNamesOfCoordSystemsRetrieved] = qKEN(c, szNamesOfCoordSystems, bufferSize)
% PI MATLAB Class Library Version 1.1.1
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.


%% Parameters

functionName = [ c.libalias, '_', mfilename];


%% Check for errors

if(c.ID < 0)
    error('The controller is not connected');
end;

if(~any( strcmp(functionName, c.dllfunctions )))
    error('%s not found',functionName);
end


%% Call the libary function

% Create temporarily variables for the C function call
bufferSize = 1000;
charBuffer  = blanks(bufferSize+1);

% helper variable to handle buffer overflow of charBuffer
bufferOverflow = true;

% Call library
while (bufferOverflow)
    try
        [bRet, ~, charBuffer] = calllib(c.libalias, functionName, c.ID, szNamesOfCoordSystems, charBuffer, bufferSize);
        if(bRet==0)
            iError = GetError(c);
            
            % Error == -5 (Buffer overflow) --> internal error handling
            if (iError == -5)
                bufsize     = bufsize * 10;
                charBuffer  = blanks(bufsize+1);
                
                if (bufsize > 100000000)
                    error('The data to receive is larger than 100 MB.');
                end
                
            % Error != -5 --> show error to user
            else
                szDesc = TranslateError(c,iError);
                error(szDesc);
            end
        else
            bufferOverflow = false;
        end
    catch
        rethrow(lasterror);
    end
end

function [bRet, charBuffer] = qKLS (c, szNameOfCoordSystem, szItem1, szItem2)
% function bRet = qKLS (c, szNameOfCoordSystem, szAxes, ValueArray)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

%% Parameters

functionName = [ c.libalias, '_', mfilename];


%% Check for errors

if (c.ID < 0)
    error ('The controller is not connected');
end

if (~any ( strcmp (functionName, c.dllfunctions )))
    error ('%s not found',functionName);
end


%% Call the libary function

% Create temporarily variables for the C function call
bufsize     = 10000;
charBuffer  = blanks (bufsize+1);

% helper variable to handle buffer overflow of charBuffer
bufferOverflow = true;

% Call library
while (bufferOverflow)
    try
        [bRet, ~, ~, ~, charBuffer] = calllib (c.libalias, functionName, c.ID, szNameOfCoordSystem, szItem1, szItem2, charBuffer, bufsize);
        if (bRet==0)
            iError = GetError (c);
            
            % Error == -5 (Buffer overflow) --> internal error handling
            if (iError == -5)
                bufsize     = bufsize * 10;
                charBuffer  = blanks (bufsize+1);
                
                if (100000000 > bufsize)
                    error ('The data to receive is larger than 100 MB.');
                end
                
            % Error != -5 --> show error to user
            else
                szDesc = TranslateError (c,iError);
                error (szDesc);
            end
        else
            bufferOverflow = false;
        end
    catch
        rethrow (lasterror);
    end
end



function bRet = CCL(c, intComandLevel, szPassWord)
% function [bRet, charBuffer] = CCL(c, intComandLevel, szPassWord)
%PI MATLAB Class Library Version 1.2.0
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

% Call library
    try
        bRet = calllib(c.libalias, functionName, c.ID, intComandLevel, szPassWord);
        if(bRet==0)
            iError = GetError(c);
            szDesc = TranslateError(c,iError);
            error(szDesc);
        end
    catch
        rethrow(lasterror);
    end
    
    

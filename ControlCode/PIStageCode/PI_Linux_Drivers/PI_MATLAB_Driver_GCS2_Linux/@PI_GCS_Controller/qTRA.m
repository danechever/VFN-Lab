function [bRet, pdValueArray] = qTRA (c, szAxes, Components)
% function bRet = qTRA(c, szNameOfCoordSystem, szAxes, ValueArray)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

%% Parameters

functionName = [ c.libalias, '_', mfilename];


%% Check for errors

if (0 > c.ID)
    error ('The controller is not connected');
end;

if (~any( strcmp( functionName, c.dllfunctions ) ) )
    error ('%s not found', functionName);
end

axesSplitted = regexp(szAxes,'[\w]+','match');
numberOfAxes = length(axesSplitted)

if (numberOfAxes  ~=  length (Components) )
    error('The number of axes and components must be the same');
end


%% Call the libary function

% Create temporarily variables for the C function call
pdComponents = libpointer ('doublePtr', Components);

ValueArray      = zeros (1, numberOfAxes );
pdValueArray    = libpointer ('doublePtr', ValueArray);


% Call library
try
    [bRet, ~, ~, pdValueArray] = calllib (c.libalias, functionName, c.ID, szAxes, pdComponents, pdValueArray);
    if (bRet==0)
        iError = GetError (c);
        szDesc = TranslateError(c, iError);
        error (szDesc);
    end
catch
    rethrow (lasterror);
end



function c = InitializeController(c)

% This code is provided by Physik Instrumente(PI) GmbH&Co.KG
% You may alter it corresponding to your needs
% Comments and Corrections are very welcome
% Please contact us by mailing to support-software@pi.ws

try
    axesIdentifierArray = qSAI_ALL(c);
catch
    axesIdentifierArray = qSAI(c);
end

c.NumberOfAxes = length(axesIdentifierArray);
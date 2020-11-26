function [dOutValues1] = qTAV(c,iInValues1)
%function [dOutValues1] = qTAV(c)
%PI MATLAB Class Library Version 1.2.0
% This code is provided by Physik Instrumente(PI) GmbH&Co.KG.
% You may alter it corresponding to your needs.
% Comments and Corrections are very welcome.
% Please contact us by mailing to support-software@pi.ws. Thank you.

if(c.ID<0), error('The controller is not connected'),end;
functionName = [ c.libalias, '_', mfilename];
if(any(strcmp(functionName,c.dllfunctions)))
	if(nargin==1)
        try
            NumberOfAnalogInputs = qTNR(c);
        catch
            NumberOfAnalogInputs = 0;
        end
		if(NumberOfAnalogInputs==0),error('Command does not work for your system. Please use GcsCommandset and GcsGetAnswer to retrieve the information'),end;
		iInValues1 = 1:min(c.NumberOfAnalogInputs,c.NumberOfAxes);
	end
	piInValues1 = libpointer('int32Ptr',iInValues1);
	nValues = length(iInValues1);
	dOutValues1 = zeros(size(iInValues1));
	pdOutValues1 = libpointer('doublePtr',dOutValues1);
	try
		[bRet,iOutValues1,dOutValues1] = calllib(c.libalias,functionName,c.ID,piInValues1,pdOutValues1,nValues);
		if(bRet==0)
			iError = GetError(c);
			szDesc = TranslateError(c,iError);
			error(szDesc);
		end
	catch
		rethrow(lasterror);
	end
else
	error('%s not found',functionName);
end

function [dOutValues1] = qGWD(c,iStart, iNumber, iTables)

% This code is provided by Physik Instrumente(PI) GmbH&Co.KG
% You may alter it corresponding to your needs
% Comments and Corrections are very welcome
% Please contact us by mailing to support-software@pi.ws

functionName = [ c.libalias, '_', mfilename];
if(strmatch(functionName,c.dllfunctions))
    piTables = libpointer('int32Ptr',iTables);
    nTables = length(iTables);
    hlen = 1000;
    header = blanks(hlen+1);
    
    ppdData = libpointer('doublePtr');
    dOutValues1 = 0;
    try
        [bRet,iTables,ppdData,header] = calllib(c.libalias,functionName,c.ID,piTables,nTables,iStart,iNumber,ppdData,header,hlen);
    catch
        rethrow(lasterror);
    end
    if(bRet==0)
        iError = GetError(c);
        szDesc = TranslateError(c,iError);
        error(szDesc);
    end
    
    headerlines = regexp(header,'\n','split');
    sampletime = -1;
    timeind = -1;
    for n = 1:length(headerlines)
        currentline = headerlines{n};
        if(~isempty(strfind(currentline,'DIM'))),nTables = str2num(currentline(strfind(currentline,'=')+1:end));end;
        if(~isempty(strfind(currentline,'NDATA'))),iNumber = str2num(currentline(strfind(currentline,'=')+1:end));end;
        if(~isempty(strfind(currentline,'SAMPLE_RATE'))),sampletime = str2double(currentline(strfind(currentline,'=')+1:end));end;
        if(~isempty(strfind(currentline,'TIME')) && ~isempty(strfind(currentline,'NAME')))
            number = currentline( (strfind(currentline,'NAME') + 4):( strfind(currentline,'=') -1));
            timeind = str2num(number) + 1;
        end
    end
    i = 0;
    
    while(i<(nTables*iNumber))
        pause(0.1);
        
        i =  GetAsyncBufferIndex(c);
        
    end
    setdatatype(ppdData,'doublePtr',nTables,iNumber);
    dOutValues1 = ppdData.Value;
    
else
    error('%s not found',functionName);
end

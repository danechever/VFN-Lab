function [dOutValues1, header] = qDDL(c,iTables,iStart,iNumber)

functionName = ['PI_', mfilename];

% Method available?
if ( ~any ( strcmp ( functionName, c.dllfunctions ) ) ), error('Method %s not found',functionName), end;


% Create variables for C interface
piTables = libpointer('int32Ptr',iTables);
nTables = 1;
ppdData = libpointer('doublePtr');
headerLength = 1000;
header = blanks(headerLength+1);

try
    [bRet,~,ppdData,header] = calllib(c.libalias, functionName, c.ID, ...
        piTables, nTables, iStart, iNumber, ppdData, header, headerLength);

    if ( 0 == bRet )
        iError = GetError ( c );
        szDesc = TranslateError ( c, iError );
        error ( szDesc );
    end
    
catch ME
    error ( ME.message );
end


headerlines = regexp(header,'\n','split');
for n = 1:length(headerlines)
    if(~isempty(strfind(headerlines{n},'DIM'))),nTables = str2num(headerlines{n}(strfind(headerlines{n},'=')+1:end));end;
    if(~isempty(strfind(headerlines{n},'NDATA'))),iNumber = str2num(headerlines{n}(strfind(headerlines{n},'=')+1:end));end;
    if(~isempty(strfind(headerlines{n},'SAMPLE_TIME'))),sampletime = str2double(headerlines{n}(strfind(headerlines{n},'=')+1:end));end;
end
i = 0;

while(i<(nTables*iNumber))
    pause(0.1);
    i =  GetAsyncBufferIndex(c);
end
setdatatype(ppdData,'doublePtr',nTables,iNumber);
dOutValues1 = ppdData.Value';
dOutValues1 = [([0:iNumber-1]*sampletime)', dOutValues1];


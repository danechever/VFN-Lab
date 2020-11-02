function c = CloseConnection(c)

functionName = [ c.libalias, '_', mfilename];
if(strmatch(functionName,c.dllfunctions))
    calllib(c.libalias,functionName,c.ID);
    c = SetDefaults(c);
else
    error(sprintf('%s not found',functionName));
end
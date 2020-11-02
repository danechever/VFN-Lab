function AssertValueVectorSize( c, Identifiers, Values )
if GetNrAxesInString(c, Identifiers) > length(Values)
    error('values vector size does not fit identifier string size');
end


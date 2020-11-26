function number = GetNrAxesInString( c, axesstring )
    lines = regexp(axesstring,'[\w-]+','match');
    number = length(lines);

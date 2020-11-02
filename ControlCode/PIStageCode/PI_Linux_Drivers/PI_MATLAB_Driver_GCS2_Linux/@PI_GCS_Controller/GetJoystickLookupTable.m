function [dOutValues] = GetJoystickLookupTable(c)
while(GcsGetAnswerSize(c))
    line = GcsGetAnswer(c);
end
dOutValues = zeros(256,2);
GcsCommandset(c,'JLT? 1 256');
while(GcsGetAnswerSize(c))
    line = GcsGetAnswer(c);
    if(strfind(line,'# END_HEADER'))
        break;
    end;
end
n = 1;
while(GcsGetAnswerSize(c))
    line = GcsGetAnswer(c);
    dOutValues(n,:) = sscanf(line, '%f')';
    n = n + 1;
end
function str = CellArrayToStringSpaceSeperated(cellArray)
% Converts Cell Array of axis names to string of spaced axes:
% Cell Array:  '1'    '2'    'X'    'Y'    'u'
% string:      1 2 X Y u
%
%
% Test behaviour with the following function calls
% clc;
% ca = {}
% returnedString = cellArrayToStringSpaceSeperated(ca)
% ca = {''}
% returnedString = cellArrayToStringSpaceSeperated(ca)
% ca = {'1'}
% returnedString = cellArrayToStringSpaceSeperated(ca)
% ca = {'1', '2', 'X', 'Y', 'u'}
% returnedString = cellArrayToStringSpaceSeperated(ca)


if (isempty(cellArray))
    str = '';
    
elseif (~iscellstr (cellArray))
    error('Wrong input datatype. Input must be a cell array of strings.');
    
else
    str = cellArray{1};
    
    for idx = 2: length(cellArray)
        str = [str, ' ', cellArray{idx}];
    end
end
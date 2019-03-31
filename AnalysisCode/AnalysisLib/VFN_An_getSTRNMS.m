function STRNMS = VFN_An_getSTRNMS(nmWpath)
% VFN_An_getSTRNMS Return all filenames, with paths, matching the
%       description within the a directory.
%
%   - The directory should be included within the provided argument. 
%   - This function basically finds all files matching the full-path name
%   - The default dir() function is at the core of this function
%   
%   STRNMS = VFN_An_getSTRNMS(nmWpath)
%     Extract the full-path filenames for all files containing 'nmWpath'
%     within the included directory. 
%     - 'nmWpath' Filenames to look for, including paths. Wildcards such as
%               '*' are allowed to find mutliple matching filenames.
%                 This should be a cell array with all filenames
%
%     Returns
%     - 'STRNMS' string array containing full path and name for all files
%               matching the specified description.
%
%   Examples:
%      STRNMS = VFN_An_getSTRNMS({'foo\*.fits'})
%         returns the filenames for all fits files within the foo directory
%      STRNMS = VFN_An_getSTRNMS({'foo\*.fits', 'bar\dat*.fits'})
%         returns the filenames for all fits within the foo directory as
%         well as all fits files within the bar\ directory that have dat in
%         their name.
%  

% Read all filenames
fls = [];
for j = 1:length(nmWpath)
    fls = [fls; dir(char(nmWpath(j)))];
end

% Create vector of string filenames
STRNMS = strings(size(fls));
for i = 1:length(STRNMS)
    STRNMS(i) = strcat(string({fls(i).folder}), filesep, string({fls(i).name}));
end
end
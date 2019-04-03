function [xax, yax, lambdaOverD] = VFN_An_simGetAxes(nm2find, indVect, an_params)
% VFN_An_simGetAxes Return axes for plotting a simulated frame in Lambda/D
%
%   - This will return an x and y vector for the axes along with the 
%           lambda/D scaling factor (in samples/(L/D))
%   - These axes are centered at the coordinates defined by indVect:
%           indVect = (row, col) = [yax, xax]
%   - Note: this function relies on VFN_An_kwdsLoad to load the keywords
%   - Note: this function reliex on VFN_An_getKwd to read a specific kwd
%
%   [xax, yax, lambdaOverD] = VFN_An_simGetAxes(nm2find, indVect, an_params)
%     Read the keywords for the file identified by nm2find and return the 
%           axes for plotting in lambda/D.
%     - 'nm2find' is the string to look for within the list. This can be a
%               partial string; does not need to be the full name. It does,
%               however, need to be a unique identifier.
%     - 'indVect' is the coordinate around which to center the axes
%               indVect = (row, col) = [yax, xax,...]
%     - 'an_params' struct containing analysis parameters. Pertinent struct
%               element for this code is an_params.STRNMS:
%               - Vector of string filenames within which to look for
%                 nm2find. This should include the full path to the file.
%     Returns
%     - '[xax, yax]' axes for plotting
%     - 'lambdaOverD' scaling factor in samples per lambda/D
%
%   Examples:
%      [xax, yax, LD] = VFN_An_simGetAxes('foo\bar', [20,22,...], an_params)
%         Returns the x,y axes and L/D value for the file containing 
%         'foo\bar' in its name. Axes will be centered at [20,22].
%

% Use values from keywords in fits header
kwds    = VFN_An_kwdsLoad(nm2find, 1, an_params);
apRad   = VFN_An_getKwd(kwds, 'aprad');
nax1    = VFN_An_getKwd(kwds, 'naxis1'); 
nax2    = VFN_An_getKwd(kwds, 'naxis2');
% Use equation from simulation code to lambda/D (L/D = sample/aprad/2)
lambdaOverD = nax1/apRad/2;     %<-- assume both axes have same L/D value
% Get axes
xax     = ((1:nax1)-indVect(2))/lambdaOverD;
yax     = ((1:nax2)-indVect(1))/lambdaOverD;

end
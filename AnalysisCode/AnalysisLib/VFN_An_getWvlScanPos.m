function [BWs, Wvls] = VFN_An_getWvlScanPos(nm2find, an_params)
% VFN_An_getWvlScanPos Return the values for the bandwidths and wavelength centers
% scanned
%   
%   - This function will calculate the wavelength and bandwidth positions
%   from a wavelength scan
%   - It will return a BWs and Wvls vectors with the positions 
%   - Note: this function relies on VFN_An_kwdsLoad to load the keywords
%   - Note: this function reliex on VFN_An_getKwd to read a specific kwd
%   
%   [BWs, Wvls] = VFN_An_getWvlScanPos(nm2find, an_params)
%     Read the keywords for the file identified by nm2find and return the 
%           positions scanned for a wavelength scan
%     - 'nm2find' is the string name to look for within the list of names. 
%               This can be a partial string; does not need to be the full 
%               name. It does, however, need to be a unique identifier.
%     - 'an_params' struct containing analysis parameters. Pertinent struct
%               elements for this code are:
%               - an_params.STRNMS - Vector of string filenames within which 
%                 to look for nm2find. This includes the full path to file.
%
%     Returns
%     - '[BWs, Wvls]' bandwidths and wavlengths scanned
%
%   Examples:
%      [BWs, Wvls] = VFN_An_getWvlScanPos('foo\bar', [20,22,...], an_params);
%         Returns the x,y axes for the cube with a filename containing 
%         'foo\bar' in its name. Axes will be centered at [20,22]

%--Finding positions
    
    %-- Load keywords
    kwds    = VFN_An_kwdsLoad(nm2find, 1, an_params);
    
    %-- Read Pertinent Keywords
    wvlmin  = VFN_An_getKwd(kwds, 'WVLMIN');
    wvlmax     = VFN_An_getKwd(kwds, 'WVLMAX');
    wvlpts     = VFN_An_getKwd(kwds, 'WVLPTS');
    bwmin  = VFN_An_getKwd(kwds, 'BWMIN');
    bwmax     = VFN_An_getKwd(kwds, 'BWMAX');
    bwpts     = VFN_An_getKwd(kwds, 'BWPTS');
    
    %-- Create vectors with every BW and Wvl value scanned
    BWs = linspace(bwmin,bwmax,bwpts);
    Wvls = linspace(wvlmin,wvlmax,wvlpts);
    
end
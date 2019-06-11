function [lower,upper] = VFN_NKT_getWvlRange(NKT)
%[lower,upper] = VFN_NKT_getWvlRange(NKT) Get wavelength range of the NKT Varia.
%
%   Inputs:
%       'NKT' - object containing NKT-related stuff
%
%   Outputs:
%       'lower' - lower wavelength (nm)
%       'upper' - upper wavelength (nm)
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019

    % Get the current lower and upper wavelengths  
    lower = NKT.nktobj.get_varia_lwpsetpoint();
    upper = NKT.nktobj.get_varia_swpsetpoint();

end

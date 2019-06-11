function VFN_NKT_setPowerLevel(NKT,powerlevel)
%VFN_NKT_setPowerLevel(NKT,powerlevel) Set power of the NKT SuperK.
%
%   Inputs:
%       'NKT' - object containing NKT-related stuff
%       'powerlevel' - power level setting (0-100)
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019

    NKT.nktobj.set_powerlevel(powerlevel);
    VFN_NKT_getPowerLevel(NKT);

end

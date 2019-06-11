function lvl = VFN_NKT_getPowerLevel(NKT)
%lvl = VFN_NKT_getPowerLevel(NKT) Get power level of the NKT SuperK.
%
%   Inputs:
%       'NKT' - object containing NKT-related stuff
%
%   Outputs:
%       'lvl' - power level setting (0-100)
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019

    lvl = NKT.nktobj.get_powerlevel();
    
    disp(['NKT emission = ',num2str(lvl)]);

end

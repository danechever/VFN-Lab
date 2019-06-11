function VFN_NKT_setEmission(NKT,on_off)
%VFN_NKT_setEmission(NKT,on_off) Toggle emission of the NKT SuperK.
%
%   Inputs:
%       'NKT' - object containing NKT-related stuff
%       'on_off' -  true: sets emission to on
%                   false: sets emission to off
%
%   author: G. Ruane
%   last modified: March 22, 2019
%   last modified: D. Echeverri June 11, 2019

    emsn = NKT.nktobj.get_emission();

    trial = 1;
    % The NKT sometimes needs to be pinged a few times before it responds
    while(and(on_off ~= emsn,trial<=NKT.numTries))
        ret = NKT.nktobj.set_emission(on_off);
        emsn = NKT.nktobj.get_emission();
        pause(NKT.delay);
        trial = trial + 1;
    end
    
    stringoptions = {'off','on'};
    disp(['     NKT emission = ',stringoptions{emsn+1}]);

    % % Save backup bench object
    % hcst_backUpBench(bench)

end

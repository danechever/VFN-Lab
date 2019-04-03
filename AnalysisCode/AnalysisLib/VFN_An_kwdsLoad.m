function reskwds = VFN_An_kwdsLoad(nm2find, mat2pull, an_params)
% VFN_An_kwdsLoad Return all keywords from the specific fits file
%
%   - This will extract all the keywords from one of the 4 VFN data cubes 
%   - NOTE: This can be used to load kwds from any fits file within 
%       an_params.STRNMS, not just VFN data cubes. Just make sure to 
%       provide a unique nm2find and SET MAT2PULL = 1.
%   
%   reskwds = VFN_An_kwdsLoad(nm2find, mat2pull, an_params)
%     Load the keywords from the file chosen by mat2pull containing the 
%     string name, nm2find,from the list of string names, an_params.STRNMS.
%     - 'nm2find' is the string to look for within the list. This can be a
%               partial string; does not need to be the full name. It does,
%               however, need to be a unique identifier.
%     - 'mat2pull' which of the 4 data cubes to read:
%               mat2pull = 2: GAINS
%               mat2pull = 3: NORMREDU
%               mat2pull = 4: SEMIREDU
%               else:         Raw data matrix
%     - 'an_params' struct containing analysis parameters. Pertinent struct
%               element for this code is an_params.STRNMS:
%               - Vector of string filenames within which to look for
%                 nm2find. This should include the full path to the file.
%     Returns
%     - 'reskwds' all keywords within the specified fits file
%
%   Examples:
%      reskwds = VFN_An_kwdsLoad(nm2find, 1, an_params)
%         Returns the keywords from the raw data-cube containing 'foo\bar' 
%         in its name from the list of full-path names, an_params.STRNMS.
%

%--Keywords laoding function
% mat2pull = 2: GAINS
% mat2pull = 3: NORMREDU
% mat2pull = 4: SEMIREDU
% else, Returns the raw data matrix.
    STRNMS = an_params.STRNMS;
    nm2find = strip(nm2find);
    inds  = contains(STRNMS, nm2find);
    flnms   = STRNMS(inds);
    switch mat2pull
        case 2
            % GAINS cube requested
            cubenm = '_GAINS';
        case 3
            % NORMREDU matrix requested
            cubenm = '_NORMREDU';
        case 4
            % SEMIREDU cube requested
            cubenm = '_SEMIREDU';
        otherwise
            % Raw data cube requested
            cubenm = {'_GAINS', '_NORMREDU', 'SEMIREDU'}';                
    end
    if mat2pull > 4 || mat2pull < 2
        % Raw data was requested
        flnmi = ~contains(flnms,cubenm);
    else
        % A different cube was requested
        flnmi = contains(flnms, cubenm);
    end
    resinfo    = fitsinfo(char(flnms(flnmi)));
    reskwds    = resinfo.PrimaryData.Keywords;      % all keywords
end
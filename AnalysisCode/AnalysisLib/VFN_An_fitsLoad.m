function resfits = VFN_An_fitsLoad(nm2find, mat2pull, an_params)
% VFN_An_fitsLoad Load a specific VFN data-cube matrix
%   
%   - This functino will load one of the 4 data-cubes output by the VFN
%     control code. 
%   - See the online experiment notes for details about each of these cubes
%   - This function uses fitsread() to read the data.
%   
%   resfits = VFN_An_fitsLoad(nm2find, mat2pull, an_params)
%     Load the cube chosen by mat2pull containing the string name, nm2find,
%     from the list of string names, an_params.STRNMS.
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
%     - 'resfits' cube of data returned by fitsread
%
%   Examples:
%      resfits = VFN_An_fitsLoad('foo\bar', 1, an_params);
%         Returns the raw data-cube containing 'foo\bar' in its name from
%         the list of full-path filenames, an_params.STRNMS.
%       

%-- Data loading function
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
    resfits    = fitsread(char(flnms(flnmi)));
end
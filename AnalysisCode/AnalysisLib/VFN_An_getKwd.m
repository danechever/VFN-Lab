function kwd = VFN_An_getKwd(kwds, kwd2find)
% VFN_An_getKwd Return the requested keyword value from within kwds
%   
%   - Note: kwd2find will automatically be capitalized to make sure it 
%       matches fits convention.
%   - Note: if kwd2find is not present in the kwds array, [] will be returned.
%
%   kwd = VFN_An_getKwd(kwds, kwd2find)
%     Extract the keyword from the set of keywords
%     - 'kwds' cell array containing keywords within which to find kwd2find
%     - 'kwd2find' keyword to find and return from the cell array
%
%     Returns
%     - 'kwd' keyword value 
%
%   Examples:
%      kwd = VFN_An_getKwd(kwds, 'NAXIS1')
%         Returns the 'NAXIS1' value from the keywords.
%  

%-- Keyword Read Function
    kwdind  = ismember(kwds(:,1), upper(kwd2find));
    kwd     = cell2mat(kwds(kwdind, 2));
end
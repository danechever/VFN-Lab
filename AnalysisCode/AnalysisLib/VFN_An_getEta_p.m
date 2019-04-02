function [eta_p, ind] = VFN_An_getEta_p(nrmDat)
% VFN_An_getEta_p Return the peak(planet) value and indices within a frame.
%   
%   - This function assumes the planet coupling is the max in the frame.
%   
%   [eta_p, ind] = VFN_An_getEta_p(nrmDat)
%     Find the planet coupling in a given frame.
%     - 'nrmDat' normalized eta-map within which to find planet coupling.
%
%     Returns
%     - 'eta_p' planet coupling (max within frame)
%     - 'ind' 6-element vector containing indices of planet coupling. 
%           NOTE: this returns y, x, z, ... coordinates (ie. row, col, ...)
%
%   Examples:
%      [eta_p, ind] = VFN_An_getEta_p(nrmDat);
%         Returns the planet coupling in the normalized data matrix, nrmDat 
%  

%--Eta_p finder - searches for max in the whole frame
    % Get minimum and linear index
    [eta_p, mxInd]  = max(nrmDat(:));
    % Translate to regular indices
    [I1, I2, I3, I4, I5, I6] = ind2sub(size(nrmDat), mxInd);
    % Format the index vector
    ind = [I1, I2, I3, I4, I5, I6];    
end
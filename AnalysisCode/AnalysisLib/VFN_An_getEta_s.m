function [eta_s, ind] = VFN_An_getEta_s(nrmDat, cropVal)
% VFN_An_getEta_s Return the star(null) coupling in the given frame
%   
%   - This function assumes the null is the minimum value within the
%       central region.
%   
%   [eta_s, ind] = VFN_An_getEta_s(nrmDat, cropVal)
%     Find the star (null) coupling in a given frame.
%     - 'nrmDat' normalized eta-map within which to find star coupling.
%     - 'cropVal' OPTIONAL: crop region to remove from analysis frame. 
%               1/cropVal of the frame is removed from all edges.
%               If ommitted: cropVal = 7 will be used.
%
%     Returns
%     - 'eta_s' star coupling (min within cropped frame)
%     - 'ind' 6-element vector containing indices of star coupling. 
%           NOTE: this returns y, x, z, ... coordinates (ie. row, col, ...)
%
%   Examples:
%      [eta_p, ind] = VFN_An_getEta_s(nrmDat);
%         Returns the star coupling in the normalized data matrix, nrmDat,
%         within the default 5/7th's region
%      [eta_p, ind] = VFN_An_getEta_s(nrmDat, cropVal);
%         Returns the star coupling in the normalized data matrix, nrmDat,
%         within the provided crop region.
%  

%--Eta_s finder - searches for minimum point within the central region
    % NOTE: this returns y, x, z , ...
    if nargin == 1
        % Use default crop window
        cropVal = 7;
    end
    %-- Crop matrix to central region in case the full donut is shown
    rowmin = max(floor(size(nrmDat, 1)/cropVal),1);
    colmin = max(floor(size(nrmDat, 2)/cropVal),1);
    rowmax = (floor((cropVal-1)*size(nrmDat,1)/cropVal+1));
    colmax = (floor((cropVal-1)*size(nrmDat,2)/cropVal+1));
    nrmDatCR = nrmDat(rowmin:rowmax,colmin:colmax,:,:,:,:);
    
    %-- Take the min of the cropped region along all dimensions 
    [eta_s, mnInd] = min(nrmDatCR(:)); %get min and linear index
    [I1, I2, I3, I4, I5, I6] = ind2sub(size(nrmDatCR), mnInd); %linear ind to subind
    
    %-- Translate submatrix indexes to full matrix ones
    rowMn     = I1+rowmin-1;
    colMn     = I2+colmin-1;
	
    %-- Format the index vector
    ind     = [rowMn, colMn, I3, I4, I5, I6];    
end
function figh = VFN_An_fitsDisp(fitsFrame, xax, yax, figN)
% VFN_An_fitsDisp Display an image using basic image parameters
%   
%   - This function will take in a data frame and display it on a white
%       background, square grid, tick marks outward, and gray colormap.
%   - It takes an optional figure handle to allow for replotting on same
%       figure or use of predefined figure parameters.
%   - Note: uses imagesc() function 
%   
%   figh = VFN_An_fitsDisp(fitsFrame, xax, yax, figN)
%     Displays fitsFrame with xax, yax axes on optional figure figN
%     Arguments
%     - 'fitsFrame' is the image (matrix) to display.
%     - 'xax', 'yax' are vectors representing the axes to use
%     - 'figN' OPTIONAL figure handle on which to plot
%     Returns
%     - 'figh' figure handle
%
%   Examples:
%      fighSamp = VFN_An_fitsDisp(frame1, xax1, yax1);
%         Displays frame1 with xax1 and yax1 axes.
%
%      fighSamp = VFN_An_fitsDisp(frame1, xax1, yax1, figh1);
%         Same as above except plots onto figh1 without changing its color
%         properties
%       

%--Basic display function
    if nargin == 1
        xax = [];
        yax = [];
    end
    if nargin == 4
        figh = figure(figN);
    else
        figh = figure('color', 'white');
    end
    imagesc(xax, yax, fitsFrame)
    axis image; axis xy; colormap gray; set(gca,'TickDir','out')
end
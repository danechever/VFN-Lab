function [theta, r] = VFN_An_fitsPol_OLD(fitsFrame, xCent, yCent)
%   This was a cartesian-to-polar coordinate conversion function I wrote.
%   It has since been replaced by the PolarTransform function from the
%   VFN-Simulations library.
%
%--Cartesian to Polar Function
    % Returns the polar coordinates for the elements in a given fitsFrame
    [xL, yL] = size(fitsFrame);
    if nargin == 3
        % center coordinates were provided; calculate cartesian shift
        shftX   = round(xL/2) - xCent; 
        shftY   = round(yL/2) - yCent;
    else
        % center coors not provided; shifts are 0
        shftX   = 0; shftY = 0;
    end
    % Create vectors for meshgrid with origin centered by shift
    X = -1+(shftX*2/(xL-1)):2/(xL-1):1+(shftX*2/(xL-1)); 
    Y = -1+(shftY*2/(yL-1)):2/(yL-1):1+(shftY*2/(yL-1));
    [x, y]  = meshgrid(X,Y);
    % Reshape into vectors
    x = x(:); y = y(:);
    [theta, r] = cart2pol(x,y);
    % It seems there are precision issues at far-off decimal places which 
        % make r, theta not symmetric about origin. Rounding solves this:
    r = round(r, 8); theta = round(theta, 8);
end
function [xax, yax] = VFN_An_getAxes(nm2find, indVect, an_params)
% VFN_An_getAxes Return the axes for plotting a frame in Lambda/D
%   
%   - This function will calculate the axes for a data frame in Lambda/D.
%   - It will return an x and y vector for the axes. 
%   - These axes are centered at the coordinates defined by indVect:
%           indVect = (row, col) = [yax, xax]
%   - Note: this function relies on VFN_An_kwdsLoad to load the keywords
%   - Note: this function reliex on VFN_An_getKwd to read a specific kwd
%   
%   nrmDat = VFN_An_getAxes(nm2find, indVect, an_params)
%     Read the keywords for the file identified by nm2find and return the 
%           axes for plotting in lambda/D.
%     - 'nm2find' is the string name to look for within the list of names. 
%               This can be a partial string; does not need to be the full 
%               name. It does, however, need to be a unique identifier.
%     - 'indVect' is the coordinate around which to center the axes
%               indVect = (row, col) = [yax, xax,...]
%     - 'an_params' struct containing analysis parameters. Pertinent struct
%               elements for this code are:
%               - an_params.STRNMS - Vector of string filenames within which 
%                 to look for nm2find. This includes the full path to file.
%               - an_params.um2LD - microns to Lambda/D conversion 
%               - an_params.V2um - volts to microns conversion (if needed)
%
%     Returns
%     - '[xax, yax]' axes for plotting (in units of lambda/D)
%
%   Examples:
%      [xax, yax] = VFN_An_getAxes('foo\bar', [20,22,...], an_params);
%         Returns the x,y axes for the cube with a filename containing 
%         'foo\bar' in its name. Axes will be centered at [20,22]

%--Axes finder - axes in L/D relative to null point
    %-- Extract pertinent parameters
    um2LD   = an_params.um2LD;
    if isfield(an_params,'V2um'); V2um = an_params.V2um; end
    
    %-- Load keywords
    kwds    = VFN_An_kwdsLoad(nm2find, 1, an_params);
    
    %-- Read Pertinent Keywords
    stepSz  = VFN_An_getKwd(kwds, 'XYSTEPS');
    ax1     = VFN_An_getKwd(kwds, 'NAXIS1');
    ax2     = VFN_An_getKwd(kwds, 'NAXIS2');
    scantype = VFN_An_getKwd(kwds, 'CNTRLCD');
    
    %-- Create axis w.r.t provided origin (indVect)
    % Get distance from origin [in native units for control code]
    xax = ((1:ax1)-indVect(2)) * stepSz; 
    yax = ((1:ax2)-indVect(1)) * stepSz; 
    % Convert to microns
    if strcmp(scantype, 'Main_Trans_PIStg')
        % Native code axis units are mm 
        % Convert axis from mm to um 
        xax = xax * 1e3; yax = yax * 1e3;
    else
        % Assume these cubes use the Thorlabs Piezo actuators so native
        % code units should be Volts
        % Convert axis from V to um
        xax = xax * V2um; yax = yax * V2um;
    end
    % Convert to L/D
    xax = xax * um2LD; yax = yax * um2LD;    
end
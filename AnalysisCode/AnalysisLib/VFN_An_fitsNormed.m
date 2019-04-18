function nrmDat = VFN_An_fitsNormed(nm2find, frame, an_params)
% VFN_An_fitsNormed Return the normalized eta-map for a scan
%   
%   - This function will load and analyze a single frame from a scan to 
%       return the normalized coupling map.
%   - The normalization includes:
%       - subtract photodiode biases (obtained from gains matrix)
%       - average the samples at each point (nominally 100 samples)
%       - divide by power incident on the fiber (measured by the red PM 
%           before the last lens and accounting for losses in that lens)
%       - account for fiber Fresnel reflections and propogation losses
%   - Note: this function relies on VFN_An_fitsLoad to load the data.
%   - Note: this function relies on VFN_An_kwdsLoad to load the keywords
%
%   * UPDATE [04/17/19]:
%       Function was updated to reflect the fact that the red PM normalization
%       value is now optional in the scan code. It now checks if the red PM was
%       used. If it was not, it does not normalize by the power on the fiber tip
%       but it does still account for propagation and fresnel losses.
%   
%   nrmDat = VFN_An_fitsNormed(nm2find, frame, an_params)
%     Load the cube/scan containing the string name, nm2find, and analyze
%       the data to return the properly normalized data.
%     - 'nm2find' is the string name to look for within the list of names. 
%               This can be a partial string; does not need to be the full 
%               name. It does, however, need to be a unique identifier.
%     - 'frame' is the frame within a scan that should be loaded. If the
%               focus and vortex were not scanned, ie. the cube is a
%               single-frame scan, frame should be [1,1,1]
%               For multi-frame scans, frame = [4th,5th,6th] dimension
%               indices in the cube. This should be foc, vortX, vortY
%               indices with the current cube structure.
%     - 'an_params' struct containing analysis parameters. Pertinent struct
%               elements for this code are:
%               - an_params.STRNMS - Vector of string filenames within which 
%                 to look for nm2find. This includes the full path to file.
%               - an_params.gnFactor - Gain factor is the scaling factor
%                 between gain settings on the photodiode.
%               - an_params.lensTHPT - fractional throughput of final lens
%               - an_params.rdOfmto - conversion from fmto to red PM
%               - an_params.FR - fractional thpt after Fresnel reflections
%               - an_params.PL - fractional thpt after propagation losses
%
%     Returns
%     - 'nrmDat' normalized eta-map as a 2D matrix
%
%   Examples:
%      nrmDat = VFN_An_fitsNormed('foo\bar', [1,1,1], an_params);
%         Returns the normalized eta-map for the single-frame cube with 
%         a filename containing 'foo\bar' in its name. 
%
%      nrmDat = VFN_An_fitsNormed('foo\bar', [3,2,4], an_params);
%         Returns the normalized eta-map for the [3,2,4] frame in the cube 
%         with a filename containing 'foo\bar' in its name. [3,2,4] are the
%         indices of the [4th,5th,6th] dimensions in the cube.
%  

%--Eta Map Calculator Function - Single Frame Scans
    %-- Define parameters
    gnFactor    = an_params.gnFactor;
    lensTHPT    = an_params.lensTHPT;
    rdOfmto     = an_params.rdOfmto;
    FR          = an_params.FR;
    PL          = an_params.PL;
    
    %-- Load the cubes
    raw     = VFN_An_fitsLoad(nm2find, 1, an_params);
    scl     = VFN_An_fitsLoad(nm2find, 2, an_params);
    kwds    = VFN_An_kwdsLoad(nm2find, 1, an_params);
    
    %-- Extract the given frame
    raw     = raw(:,:,:,frame(1),frame(2),frame(3));
    scl     = scl(:,:,:,frame(1),frame(2),frame(3));
    
    %-- Check vort/foc were not scanned - ie. effectively a single scan
    if VFN_An_getKwd(kwds, 'NAXIS') ~= 3
        error('Matrix dimensionality is not correct for this analysis:\n    %s', nm2find)
    end
    
    %-- Subtract bias from raw data
    % We know these data matrices only have 3 dimensions so just use repmat
    gns     = repmat(scl(:,:,3),[1 1 size(raw,3)]);    %bias matrix
    % Subtract element-wise now that dimensionality is matched
    nrm     = raw - gns;
    
    %-- Calculate the mean of the Nread samples for
    nrm = mean(nrm,3);
    
    %-- Scale the data back to femto setting 10^6
    gns = gnFactor*ones(size(scl(:,:,1)));  % Gain matrix
    gns = gns.^-(scl(:,:,1)-6);
    nrm = nrm.*gns;
    
    %-- Check if red PM was used. Assumes that if normflag keyword is missing,
        %the cube predates this keyword and has the red PM data in it.
        % Note: getKwd() returns [] if the keyword was not present
    pmflg = VFN_An_getKwd(kwds, 'NORMFLAG');
    if isempty(pmflg) || pmflg == 'T'
        %-- Get red thorlabs PM value for normalization since it exists
        % Read from keywords
        nrmVl = VFN_An_getKwd(kwds, 'NORMMEAN');
        % Account for losses in the final lens
        nrmVl = nrmVl*lensTHPT;
        % Convert from Watts to uWatts
        nrmVl = nrmVl*10^6;
    else
        %-- Red PM was not used (flag was present but not 'T[rue]'
        % Set nrmVal to 1 so that we do not scale/normalize the data
        nrmVl = 1;
    end
    
    %-- Normalize the data to get eta's (use nem's equation)
    nrmDat = nrm.*(rdOfmto*(1/nrmVl)*(1/(FR^2))*(1/(PL^2)));
end
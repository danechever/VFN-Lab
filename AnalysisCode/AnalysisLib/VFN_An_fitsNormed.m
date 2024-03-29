function nrmDat = VFN_An_fitsNormed(nm2find, frame, an_params)
% VFN_An_fitsNormed Return the normalized eta-map for a scan
%   
%   - This function will load and analyze a single frame from a scan to 
%       return the normalized coupling map.
%   - The normalization includes:
%       - subtract photodiode biases (obtained from gains matrix)
%       - average the samples at each point (nominally 100 samples)
%       - divide by power incident on the fiber
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
%   * UPDATE [07/2021]:
%       Correct for when the redPM is behind the lens s.t. lens loss is not
%       accounted for.
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
    FR          = an_params.FR;
    PL          = an_params.PL;
    % Make lens throughput optional since redPM is sometimes behind lens
    if isfield(an_params,'lensTHPT'); lensTHPT = an_params.lensTHPT; end
    
    %-- Load the cubes
    raw     = VFN_An_fitsLoad(nm2find, 1, an_params);
    scl     = VFN_An_fitsLoad(nm2find, 2, an_params);
    kwds    = VFN_An_kwdsLoad(nm2find, 1, an_params);
    
    %-- Extract the given frame
    raw     = raw(:,:,:,frame(1),frame(2),frame(3));
    scl     = scl(:,:,:,frame(1),frame(2),frame(3));
    
    % Extract the type of scan
    scantype = VFN_An_getKwd(kwds, 'CNTRLCD');
    
    %-- For wavelength scans, check that vort/foc were not scanned
      % ie. effectively a single scan
    if ((VFN_An_getKwd(kwds, 'NAXIS') ~= 3) && (scantype ~= 'WavelengthScan'))
        error('Matrix dimensionality is not correct for this analysis:\n    %s', nm2find)
    end
    
    %-- Subtract bias (frame 3 in scl matrix) from raw data
    % *matlab automatically resizes compatible matrices for subtraction*
    nrm     = raw - scl(:,:,3);
    
    %-- Calculate the mean of the Nread samples for
    nrm = mean(nrm,3);
    
    %-- Scale the data back to femto setting 10^6
    gns = gnFactor*ones(size(scl(:,:,1)));  % Gain matrix
    gns = gns.^-(scl(:,:,1)-6);
    nrm = nrm.*gns;
    
    
    if strcmp(scantype, 'WavelengthScan')
        rdOfmto     = VFN_getFmtoResponsivity(an_params.Wvls(frame(1)));%an_params.rdOfmto;
        if isempty(an_params.NormVals)
            nrmvals = nan(length(an_params.Wvls), length(an_params.BWs));
            for j = 1:length(an_params.Wvls)
                for i = 1:length(an_params.BWs)
                    nrmvals (j,i) = VFN_An_getKwd(kwds, ['RAWNRM' sprintf('%02d', (j-1)*length(an_params.BWs)+i)]);
                end
            end
            if isempty(nrmvals(1,1))
                for j = 1:length(an_params.Wvls)
                    for i = 1:length(an_params.BWs)
                        nrmvals (j,i) = VFN_An_getKwd(kwds, ['NRMVAL' sprintf('%02d', (j-1)*length(an_params.BWs)+i)]);
                    end
                end
                nrmVl = nrmvals(frame(1),frame(2));
            else
                nrmVl = nrmvals(frame(1),frame(2))/VFN_getFmtoResponsivity(an_params.Wvls(frame(1)));
            end
        else
            nrmVl = an_params.NormVals(frame(1),frame(2));
        end 
    else
    % This "else" includes 'Main_Trans_PIStg'
        rdOfmto     = an_params.rdOfmto;
        %-- Check if red PM was used. Assumes that if normflag keyword is missing,
            %the cube predates this keyword and has the red PM data in it.
            % Note: getKwd() returns [] if the keyword was not present
        pmflg = VFN_An_getKwd(kwds, 'NORMFLAG');
        if isempty(pmflg) || pmflg == 'T'
            %-- Get red thorlabs PM value for normalization since it exists
            % Read from keywords
            nrmVl = VFN_An_getKwd(kwds, 'NORMMEAN');
            %NOTE: Main_Trans_PIStg doesn't need to account for lens losses
            % since the redPM was behind the lens.
            if ~strcmp(scantype,'Main_Trans_PIStg')
                % Account for losses in the final lens
                nrmVl = nrmVl*lensTHPT;
            end
            % Convert from Watts to uWatts
            nrmVl = nrmVl*10^6;
            % Convert from uWatts to Femto volts @ gain=6
            nrmVl = nrmVl/rdOfmto;
        else
            %-- Red PM was not used (flag was present but not 'T[rue]'
            % Set nrmVal to 1 so that we do not scale/normalize the data
            nrmVl = 1;%/rdOfmto;
        end
    end
    
    %-- Normalize the data to get eta's (use nem's equation)
    nrmDat = nrm.*((1/nrmVl)*(1/(FR^2))*(1/(PL^2)));
end
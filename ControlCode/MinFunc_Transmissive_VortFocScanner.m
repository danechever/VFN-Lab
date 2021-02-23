function relTint = MinFunc_Transmissive_VortFocScanner(X, PIdevs, Zabs, PM, dist, pmNread, isPMNorm, pmCalWvl, zabBackPos)

global itr nl_min_hists FMTO

%% Prep motion
%-- Extract desired positions from vector
fibZ  = X(1);    % Fiber X position [in mm]
vortX = X(2);    % Fiber Y position [in mm]
vortY = X(3);    % Fiber Z position [in mm]

flg = false;
%-- Check that all values are within allowable range
if ~(PIdevs.fibZ_lower<=fibZ && fibZ<=PIdevs.fibZ_upper)
    % Z fiber position is beyond range
    relTint = 1000; 
    flg = true;
end
if ~(Zabs.vortX_lower<=vortX && vortX<=Zabs.vortX_upper)
    % X vortex position is beyond range
    relTint = 1000;
    flg = true;
end
if ~(Zabs.vortY_lower<=vortY && vortY<=Zabs.vortY_upper)
    % Y vortex position is beyond range
    relTint = 1000; 
    flg = true;
end

if flg
    %-- Save values in history struct
    % Increment itr counter
    itr = itr + 1;
    nl_min_hists.X(itr,:)  = [fibZ, vortX, vortY];
    nl_min_hists.PWR(itr)  = relTint;
    nl_min_hists.nullPt(itr,:) = [nan, nan];
    nl_min_hists.nullVl(itr) = nan;
    nl_min_hists.peakVl(itr) = nan;
    nl_min_hists.normVl(itr) = nan;
    return
end

%% Move
% Move in fiber Z
VFN_PIStage_move(PIdevs.fibZ, fibZ);

%-- Move in vortex X
% Move to zabBackPos
VFN_Zab_move(Zabs.vortX, zabBackPos(1));
% Move to desired pos
VFN_Zab_move(Zabs.vortX, vortX);

%-- Move in vortex Y
% Move to zabBackPos
VFN_Zab_move(Zabs.vortY, zabBackPos(2));
% Move to desired pos
VFN_Zab_move(Zabs.vortY, vortY);

%- Make sure pmX is out of the way
if isfield(Zabs, 'pmX')
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
elseif isPMNorm
    error('PM Norm is requested but no pmX Zaber exists.')
end

%-- Wait for power to settle (likely not needed but oh well)
pause(0.2);

%% Do course raster scan to get donut
[meas, scales, measScl, FMTO] = VFN_Raster_Scan_XY(dist, PIdevs, FMTO, 0.1, false);

%% Measure power at Fiber Tip (for normalization)
%-- Use the red power meter to measure the power
if isPMNorm
    %-- Move redPM into the beam
    %fprintf('\nredPM norm used. Moving fibZ...')
    
    %Move fiber stage out of the way to avoid collision with pmX
    VFN_PIStage_move(PIdevs.fibZ, PIdevs.fibZ_lower);
    %fprintf('  fibZ at %0.2f',VFN_PIStage_getPos(PIdevs.fibZ))
    
    % Brief pause just to be safe
    pause(2);
    
    %Move if no collision
    %fprintf('\nMoving pmX...')
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_upper);
    %fprintf(' Doing measurments...')

    %-- Read from red PM 
    % Wait for PM value to settle
    pause(2);

    % Take measurements of the power (use try catch to avoid errors)
    NTries = 5;
    for i = 1:NTries
        try 
            % This is where the read often fails 
            pmRead = VFN_PM100D_read(PM, pmNread);
        catch ME
            % An error occurred, let the user know
            warning('PM read failed')
            if i < 5
                % We have tried less than NTries times, try again
                
                % Wait a moment to give the redPM time to figure itself out
                pause(1);
                
                % First, try playing with the wavelength settings
                    % This seems to help for some reason
                VFN_PM100D_setWavelength(PM, pmCalWvl);
                pause(0.1);
                VFN_PM100D_getWavelength(PM);
                
                % Now try again (continue will jump straight to next iter.)
                fprintf('Trying again\n')
                continue
            else
                % We've tried 5 times, it's time to give up
                error('Multiple redPM reads failed')
            end
        end
        % If we get here, succeeded in reading so exit for-loop and proceed
        break
    end
    
    % Move redPM zaber, pmX, out of beam again
    VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
    %fprintf('\nDone. pmX moved to %0.2f',VFN_Zab_getPos(Zabs.pmX))
else
    %-- Set pmRead to -9999 if PM is not to be read
    pmRead = ones(pmNread,1)*-9999;
end


%% Normalize

%** Bias is already subtracted in the main for loop

%-- Define the propagation and fresnel reflection loss values
FL = 0.9654;  % (1-loss) Fraction per surface
PL = 0.9966;  % (1-loss) Fraction per meter of fiber

%-- Define scaled matrix; equal to mean of meas unless gain was varied
measScl = mean(measScl,3).*scales(:,:,2);

%-- Get FMTO sensitivity at given wavelength
%fmtoSnsScl = VFN_getFmtoResponsivity(pmCalWvl);
fmtoSnsScl = 0.62; % value from manual comparison

%-- Define a separate matrix which normalizes the data by power on the 
if isPMNorm
    %-- When actual PM reading was done, normalize as usual
    % Account for loss from fiber lens
        % 3.4% reflection spec at 780nm, 0.3% reflection measured at 635nm
    %pmRead1 = mean(pmRead)*(1-0.034);%0.997;    % Updated on 2/8/19 based on 10/26/18
    % Convert the calibration PM value from W to uW
    pmRead1 = mean(pmRead)*10^6;
    % Translate the calibration PM value from uW to V (at 10^6)
    pmRead1 = pmRead1*fmtoSnsScl;     % V/uW conversion is ~ 0.63 @633; ~0.9 @ 780
    % Account for Fresnel and Propagation losses
    nrmVal = (1/pmRead1)*(1/(FL^2))*(1/(PL^2));
else
    %-- When PM was not read, set normalization value to 1
    nrmVal = 1;
end

% Create the matrix and normalize by total power on fib tip
measScl2 = measScl*nrmVal;

%% Compute rel Tint

%-- Find null in central region
cropVal = 7;
[mnVal, mnInd] = VFN_An_getEta_s(measScl, cropVal);
eta_s = mnVal*nrmVal;

%-- Find eta_p as radially-averaged max 
    % measScl*nrmVal to get fully normed value directly
[radProf, ~] = VFN_An_radAverage(measScl*nrmVal,[mnInd(2), mnInd(1)]);
radMax = max(radProf);
%fprintf('Max, normed (radial average): %f\n', radMax);

% Compute Rel Tint
relTint = eta_s/(radMax^2);

%-- Save values in history struct
% Increment itr counter
itr = itr + 1;
nl_min_hists.X(itr,:)  = [fibZ, vortX, vortY];
nl_min_hists.PWR(itr)  = relTint;
nl_min_hists.nullPt(itr) = mnInd;
nl_min_hists.nullVl(itr) = eta_s;
nl_min_hists.peakVl(itr) = radMax;
nl_min_hists.normVl(itr) = nrmVal;
end
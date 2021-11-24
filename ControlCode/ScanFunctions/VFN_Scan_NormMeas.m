function pmRead = VFN_Scan_NormMeas(PIdevs,Zabs,PM,pmNread,pmCalWvl)
% VFN_Scan_NormMeas Function for measuring the power incident on the SMF
%   during a scan
%   
% This function was written considering only the requirements of the
% VFN_Main_... scripts. As such, check it carefully before using it
% anywhere else.
%
%
% NOTE::: This function only does a redPM measurement; there is no option
% to norm with the MMF.
%
%
% Value returned in [mm]
%   
%   EXAMPLE:______
%   pos = VFN_PIStage_getPos(stage)
%       'stage' = instance of PI device (USB connection) object
%                 (ie. a single field in the PIdevs struct)
%       'pos'   = output position in [mm]

%% Measure power at Fiber Tip
%-- Move redPM into the beam
fprintf('\nredPM norm used. Moving fibZ...')

%Move fiber stage out of the way to avoid collision with pmX
VFN_PIStage_move(PIdevs.fibZ, PIdevs.fibZ_lower);
fprintf('  fibZ at %0.2f',VFN_PIStage_getPos(PIdevs.fibZ))

% Brief pause just to be safe
pause(2);

%Move if no collision
fprintf('\nMoving pmX...')
VFN_Zab_move(Zabs.pmX, Zabs.pmX_upper);
fprintf(' Doing measurments...')

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

fprintf(' Read: %f [uW].', mean(pmRead)*10^6);

% Move redPM zaber, pmX, out of beam again
fprintf('\nMoving pmX out...')
VFN_Zab_move(Zabs.pmX, Zabs.pmX_lower);
fprintf('Done. pmX moved to %0.2f\n',VFN_Zab_getPos(Zabs.pmX))
end
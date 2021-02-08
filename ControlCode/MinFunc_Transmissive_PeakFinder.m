function pwr = MinFunc_Transmissive_PeakFinder(X, PIdevs)

global itr pk_min_hists FMTO

%% Prep motion
%-- Extract desired positions from vector
fibX = X(1);    % Fiber X position [in mm]
fibY = X(2);    % Fiber Y position [in mm]
fibZ = X(3);    % Fiber Z position [in mm]

flg = false;
%-- Check that all values are within allowable range
if ~(PIdevs.fibX_lower<=fibX && fibX<=PIdevs.fibX_upper)
    % X position is beyond range
    pwr = 1000; 
    flg = true;
end
if ~(PIdevs.fibY_lower<=fibY && fibY<=PIdevs.fibY_upper)
    % Y position is beyond range
    pwr = 1000;
    flg = true;
end
if ~(PIdevs.fibZ_lower<=fibZ && fibZ<=PIdevs.fibZ_upper)
    % Z position is beyond range
    pwr = 1000; 
    flg = true;
end

if flg
    %-- Save values in history struct
    % Increment itr counter
    itr = itr + 1;
    pk_min_hists.X(itr,:)  = [fibX, fibY, fibZ];
    pk_min_hists.PWR(itr)  = pwr;
    pk_min_hists.SCLS(itr) = nan;
    return
end

%% Move
% Move in X
VFN_PIStage_move(PIdevs.fibX, fibX);
% Move in Y
VFN_PIStage_move(PIdevs.fibY, fibY);
% Move in Z
VFN_PIStage_move(PIdevs.fibZ, fibZ);

%-- Wait for power to settle (likely not needed but oh well)
pause(0.2);

%% Measure power
% Take a sample reading at current power
read    = VFN_FMTO_readN(FMTO);

% Perform auto scaling if desired
if FMTO.isAutoScale
    % Save old FMTO_scale for comparison after autoGain
    old_scale = FMTO.FMTO_scale;
    % Modify gain accordingly
    FMTO = VFN_FMTO_setAutoGain(FMTO, mean(read));
    % Check if re-read is needed
    if old_scale ~= FMTO.FMTO_scale
        % FMTO_scale changed so the gain changed
        fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO.FMTO_scale)
        read    = VFN_FMTO_readN(FMTO);
%         if FMTO.FMTO_scale > 9
%             warning('Gain >9: %i',FMTO.FMTO_scale)
%         end
    end
end
if find(read>10)
    warning('Power is too high')
end

pause(0.1)

%% Scale and format power
scales = FMTO.gnFact^-(FMTO.FMTO_scale-6);

%-- Bias subtract
locBias = VFN_FMTO_getBias(FMTO.FMTO_scale);

%-- Take multiple samples at current location and settings
%- Preallocate power vector
pwr = nan(10,1);
for i = 1:length(pwr)
    pause(0.05);
    read    = VFN_FMTO_readN(FMTO);
    pwrtmp  = read - locBias;
    %- Average and apply gain factor
    pwr(i) = mean(pwrtmp)*scales;
end
pwr = -1*mean(pwr);

%-- Save values in history struct
% Increment itr counter
itr = itr + 1;
pk_min_hists.X(itr,:)  = [fibX, fibY, fibZ];
pk_min_hists.PWR(itr)  = -1*pwr;
pk_min_hists.SCLS(itr) = scales;
end
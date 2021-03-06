function pwr = MinFuncNullPolar_FullScan(X, MDT, cent)
% This only searches for the null in fiber X and Y positioning
% ALSO, it uses polard coordinates in X so that bounds can be circular

global s itr nl_min_hists FMTO_scale RBnd

%% Prep motion
%-- Extract desired positions from vector
fibR = X(1);                % Fiber R position
fibT = X(2);                % Fiber Theta position

%-- Convert from polar to cartesian 
fibX = fibR*cos(fibT);
fibY = fibR*sin(fibT);

%-- Translate to absolute position in axis units
fibX = fibX + cent(1);
fibY = fibY + cent(2);

flg = false;
if fibR > RBnd
    pwr = 1000; 
    flg = true;
end
%-- Check that all values are within physical range
if ~(0<=fibX && fibX<=150)
    % X-piezo position is beyond range
    pwr = 1000; 
    flg = true;
    %error('Desired fibX position (%f) is beyond allowable range', fibX)
end
if ~(0<=fibY && fibY<=150)
    % Y-piezo position is beyond range
    pwr = 1000; 
    flg = true;
    %error('Desired fibY position (%f) is beyond allowable range', fibY)
end

if flg
    %-- Save values in history struct
    % Increment itr counter
    itr(2) = itr(2) + 1;
    nl_min_hists.X(itr(1),itr(2),:)   = [fibX, fibY, fibR, fibT];
    nl_min_hists.PWR(itr(1), itr(2))  = pwr;
    nl_min_hists.SCLS(itr(1), itr(2)) = nan;
    return
end

%% Move
%-- Move the Piezos
% Remove backlash by going back to 5V every time
DE2_MDTVol(MDT, 5, 'x', 0);
DE2_MDTVol(MDT, 5, 'y', 0);
% Move in X
DE2_MDTVol(MDT, fibX, 'x', 0);
% Move in Y
DE2_MDTVol(MDT, fibY, 'y', 0);

%-- Wait for power to settle (likely not needed but oh well)
pause(0.1);

%% Measure power
%-- Measure null at new location
% Take a sample reading at current power
read    = startForeground(s);
% Save old FMTO_scale for comparison after autoGain
old_scale = FMTO_scale;
% Save old read value for comparison after autoGain
old_read = mean(read);
% Modify gain accordingly
[FMTO_scale, s] = VFN_FMTO_setAutoGain(s, old_read, FMTO_scale);
% Check if re-read is needed
if old_scale ~= FMTO_scale
    % FMTO_scale changed so the gain changed
    fprintf('\nFemto gain was changed from %i to %i\n',old_scale, FMTO_scale)
    %read    = startForeground(s);
    %ratio   = ratio * (old_read/mean(read));
    if FMTO_scale > 9
        warning('Gain >9: %i',FMTO_scale)
    end
end

%% Scale and format power
%-- GAIN Factor for femto
gnFact  = 9.97;
scales = gnFact^-(FMTO_scale-6);

%-- Bias subtract
switch FMTO_scale
    case 5
        locBias = 0.004905;%0.0033887;
    case 6
        locBias = 0.005382;%0.0042960;
    case 7
        locBias = 0.005405;%0.0044690;
    case 8
        locBias = 0.005350;%0.0047430;
    case 9
        locBias = 0.004941;%0.0068927;
    case 10
        locBias = -0.001120;
    case 11
        locBias = -0.059031;
    otherwise
        warning('No bias forFMTO_scale = %i\n',FMTO_scale)
        locBias = nan;
end

%-- Take multiple samples at current location and settings
%- Preallocate power vector
pwr = nan(10,1);
for i = 1:length(pwr)
    pause(0.05);
    read    = startForeground(s);
    pwrtmp  = read - locBias;
    %- Average and apply gain factor
    pwr(i) = mean(pwrtmp)*scales;
end
pwr = mean(pwr);

%-- Save values in history struct
% Increment itr counter
itr(2) = itr(2) + 1;
nl_min_hists.X(itr(1),itr(2),:)   = [fibX, fibY, fibR, fibT];
nl_min_hists.PWR(itr(1), itr(2))  = pwr;
nl_min_hists.SCLS(itr(1), itr(2)) = scales;
end
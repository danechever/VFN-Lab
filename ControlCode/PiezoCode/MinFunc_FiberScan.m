function pwr = MinFunc_FiberScan(X, fibZ, MDT, FMTO_scale, nrmFac)
% NOTE: backlash is always removed from zaber motion. Even when moving forward.

global s

%% Prep motion
%-- Extract desired positions from vector
fibX = X(1);                % Fiber X position
fibY = X(2);                % Fiber Y position
focZ = X(3)/nrmFac;            % Fiber Z position [IN MILLIMETERS]

%-- Check that all values are within physical range
if ~(0<=fibX && fibX<=150)
    % X-piezo position is beyond range
    error('Desired fibX position (%f) is beyond allowable range', fibX)
end
if ~(0<=fibY && fibY<=150)
    % Y-piezo position is beyond range
    error('Desired fibY position (%f) is beyond allowable range', fibY)
end
if ~(0<=focZ && focZ<=7.8)
    % Z-Zaber position is beyond range
    error('Desired focZ position (%f) is beyond allowable range', focZ)
end

%% Move
%-- Move the Piezos
% Move in X
DE2_MDTVol(MDT, fibX, 'x', 0);
% Move in Y
DE2_MDTVol(MDT, fibY, 'y', 0);

%-- Move the Zaber
% Remove backlash
VFN_Zab_move(fibZ, focZ-0.005);
% Now move
VFN_Zab_move(fibZ, focZ);

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
    read    = startForeground(s);
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
pwr = read - locBias;

%-- Average and apply gain factor
pwr = mean(pwr)*scales;
end
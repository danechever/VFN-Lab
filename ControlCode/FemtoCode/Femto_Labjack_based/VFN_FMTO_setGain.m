function VFN_FMTO_setGain(FMTO, isBlock)
% VFN_FMTO_setGain Function to change the Femto gain setting:
%
%   - Assumes the gain setting is set to "high-speed" on the Femto
%   - Includes an optional pause for the rise time at the given gain setting
%
%   Modified from HCST-R set gain function
%
%   EXAMPLE:______
%   VFN_FMTO_setAutoGain(FMTO, isBlock)
%       FMTO: a struct representing the Femto object and pertinent params     
%       isBlock: (OPTIONAL - default = true)
%                   true = pause to account for rise time at the new gain
%                   false = do not pause; return as soon as gain is set
%

%% Check if isBlock is provided
if nargin <2
    % User didn't provide isBlock explicitly so default it to true
    isBlock = true;
end

%% Check that the commanded scale setting is within bounds
if FMTO.FMTO_scale < 5 || FMTO.FMTO_scale > 11
    % beyond bounds: throw error
    error('Desired FMTO_scale is beyond bounds: %i', FMTO.FMTO_scale)
end

%% Convert gain value to binary vector so we can apply it to the DIO pins
gain_pin_dec_setting = FMTO.FMTO_scale - 5;
gain_pin_bin_setting = dec2bin(gain_pin_dec_setting, 3);
gain_pin_list = py.list();

for i=1:length(gain_pin_bin_setting)
    gain_pin_list.append(int64(str2double(gain_pin_bin_setting(i))));
end

% append high speed gain setting
gain_pin_list.append(int64(0)); 

%% Set the gain
names = py.list(FMTO.DIOPorts);
numFrames = int64(length(names));

py.labjack.ljm.eWriteNames(FMTO.LJ, numFrames, names, gain_pin_list);

%% Pause for the power to settle at the given setting
% Based on a characterization done on 6/28/2021, the photodiode+Labjack
% system takes a moment to stabilize after changing the gain. The
% settling/rise time was computed in that charac. and is applied here

if isBlock
    % ENTER THE SETTLING TIMES FOR EACH GAIN VALUE BELOW
    %gain:  5    6     7     8     9     10    11
    RT = [0.35, 0.35, 0.35, 1.00, 1.50, 1.75, 2.00];

    % Subtract 4 from FMTO_scale so we can index into the RT vector directly
    pause(RT(FMTO.FMTO_scale - 4));
end
end
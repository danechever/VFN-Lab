function VFN_FMTO_setGain(FMTO)
% VFN_FMTO_setGain Function to change the Femto gain setting:
%
%   - Assumes the gain setting is set to "high-speed" on the Femto
%
%   Modified from HCST-R set gain function
%
%   EXAMPLE:______
%   VFN_FMTO_setAutoGain(FMTO)
%       FMTO: a struct representing the Femto object and pertinent params     
%       

%% Check Validity
% TODO: check type of s with ljm library

% Check that the commanded scale setting is within bounds
if FMTO.FMTO_scale < 5 || FMTO.FMTO_scale > 11
    % beyond bounds: throw error
    error('Desired FMTO_scale is beyond bounds: %i', FMTO.FMTO_scale)
end

gain_pin_dec_setting = gainscale - 5;
gain_pin_bin_setting = dec2bin(gain_pin_dec_setting, 3);
gain_pin_list = py.list();

for i=1:length(gain_pin_bin_setting)
    gain_pin_list.append(int64(str2double(gain_pin_bin_setting(i))));
end

% append high speed gain setting
gain_pin_list.append(int64(0)); 

names = py.list({"FIO2", "FIO1", "FIO0", "FIO3"});
numFrames = int64(length(names));

py.labjack.ljm.eWriteNames(FMTO.LJ, numFrames, names, gain_pin_list);

end
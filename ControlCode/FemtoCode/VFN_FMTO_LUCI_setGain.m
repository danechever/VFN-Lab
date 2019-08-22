function code = VFN_FMTO_LUCI_setGain(FMTO_scale)
% VFN_FMTO_setGain Function to change the Femto gain setting:
%
%   - FMTO_scale must be within the bounds of the power meter:
%       [5 <= FMTO_scale <= 11]
%
%   - Assumes the gain setting is set to "high-speed" on the Femto
%   
%   - Returns the error code from the LUCI library
%
%   EXAMPLE:______
%   ERROR_CODE = VFN_FMTO_setAutoGain(FMTO_scale)
%       ERROR_CODE:          Error code when changing gain
%       FMTO_scale: Desired scale setting defined by: 10^FMTO_scale
verbose = false;
%% Check Validity
% Check that the commanded scale setting is within bounds
if FMTO_scale < 5 || FMTO_scale > 11
    % beyond bounds: throw error
    error('Desired FMTO_scale is beyond bounds: %i', FMTO_scale)
end
%% Perform change

% Translate new gain setting to binary
cmd = dec2bin(FMTO_scale-5,3);

%setting high speed mode, DC, gain 
data_low = ['00001' cmd];
data_high = '00000000';

%call dll if not already called
if not(libisloaded('LUCI_10_x64'))
    loadlibrary('LUCI_10_x64', 'LUCI_10.h')
    calllib('LUCI_10_x64', 'EnumerateUsbDevices');
end

%Show library functions
%libfunctions('LUCI_10_x64')

%testing if luci is receiving function calls; flashes led on and off
%calllib('LUCI_10_x64','LedOn',int16(1))
%pause(1)
%calllib('LUCI_10_x64','LedOff',int16(1))

code = calllib('LUCI_10_x64', 'WriteData', int32(1), int32(bin2dec(data_low)), int32(bin2dec(data_high)));

if verbose == true
    if code == 0
        disp('Femto gain changed successfully')
    elseif code == -1
        disp('Femto gain error: selected LUCI not in list')
    elseif code == -2
        disp('Femto gain error: LUCI doesnt respond')
    end
end


    
% New gain setting is done; replicate original session
%VFN_setUpFMTO
end
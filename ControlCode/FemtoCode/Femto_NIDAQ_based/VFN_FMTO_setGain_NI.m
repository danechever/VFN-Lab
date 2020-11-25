function s = VFN_FMTO_setGain(s, FMTO_scale)
% VFN_FMTO_setGain Function to change the Femto gain setting:
%
%   - FMTO_scale must be within the bounds of the power meter:
%       [5 <= FMTO_scale <= 11]
%
%   - Assumes the gain setting is set to "high-speed" on the Femto
%
%   - The daq session (s) needs to be modified and then restarted so the
%       the code tries to return a session with the same parameters.
%           *** You can add parameters at the beginning of the "Perform
%           Change" section.
%
%   EXAMPLE:______
%   s = VFN_FMTO_setAutoGain(s, FMTO_scale)
%       s:          An NI daq session 
%       FMTO_scale: Desired scale setting defined by: 10^FMTO_scale
%

%% Check Validity
% Check that s variable is correct
if ~isa(s, 'daq.ni.Session')
    error('The first argument must be an ni daq session')
end

% Check that the commanded scale setting is within bounds
if FMTO_scale < 5 || FMTO_scale > 11
    % beyond bounds: throw error
    error('Desired FMTO_scale is beyond bounds: %i', FMTO_scale)
end

%% Perform change

%-- Save the current session parameters 
%   (You can add more parameters here as needed for your application. 
%   These should be whatever you are defining when VFN_setUpFMTO gets called)
Nrate = s.Rate;             % [samples/second]
Nread = s.NumberOfScans;    % Number of samples at a given location

%Modify current session to allow output
s.Rate = 10; s.DurationInSeconds = 0.2;
% Add digitial channel (suppressing known warning)
warning('off', 'daq:Session:onDemandOnlyChannelsAdded')   
addDigitalChannel(s, 'Dev1', 'Port0/Line0:4', 'OutputOnly');
warning('on', 'daq:Session:onDemandOnlyChannelsAdded')

% Translate new gain setting to binary
cmd = de2bi(FMTO_scale-5,3);

% Change gain
outputSingleScan(s, [cmd,1,0])

% New gain setting is done; replicate original session
VFN_setUpFMTO
end
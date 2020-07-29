function FMTO = VFN_FMTO_LUCI_setAutoGain(FMTO, readVal)
% VFN_FMTO_LUCI_SETAUTOGAIN Function for automatically setting the gain:
%   Takes in a power reading and changes the gain accordingly to ensure
%   that future reads are within reasonable bounds on the Femto PM. If 
%   the reading is already within bounds, nothing is done.
%
%   - Assumes the gain setting is set to "high-speed" on the Femto
%
%   - The code shoots to intelligently place the new reading within bounds:
%       If readVal < modBoundLow, the gain will be increased until 
%         newVal > 1  
%       If readVal > modBoundHigh, the gain will be decreased until 
%         newVal < 1
%   
%       Current Mod Bounds:     [0.1 V < readVal < 9.5 V]
%
%   EXAMPLE:______
%   [FMTO] = VFN_FMTO_LUCI_SETAUTOGAIN(FMTO, readVal)
%       readVal:    A recent reading to use as a starting point 
%       FMTO:       struct of Femto parameters. Needs field:
%                   - FMTO_scale: Current scale setting defined by 10^FMTO_scale
%                   - s:          An NI daq session

modBoundLow     = 0.1; % Volts
modBoundHigh    = 9.5; % Volts

%% Determine if change is needed
% Check that FMTO.s variable is correct
if ~isa(FMTO.s, 'daq.ni.Session')
    error('The first argument must be an ni daq session')
end

if readVal > modBoundLow && readVal < modBoundHigh
    % Do nothing since within bounds
    return
end

%% Perform change

if readVal < modBoundLow
    % Initial gain was too low
    while readVal < .8
        if FMTO.FMTO_scale >= 11
            % Gain is already at maximum
            if readVal < modBoundLow
                % Warn if power is still outside of bounds
                warning(['Femto gain is at Maximum but power is still low\n' ...
                      ' Current Gain = ~10^%i   |  Current Power = %f'], ...
                      FMTO.FMTO_scale, readVal)
            end
            % Exit function
            return
        end
        
        % Increase the gain
        FMTO.FMTO_scale = FMTO.FMTO_scale + 1;
        
        % Change gain
        VFN_FMTO_LUCI_setGain(FMTO);
        
        % Check power at new gain 
        read = nan(10,1);
        for i = 1:length(read)
            read(i) = FMTO.s.inputSingleScan;
        end
        readVal = mean(read(:));        
    end
    % New gain setting is done; exit  
    return
end
if readVal > modBoundHigh
    % Initial gain was too high
    while readVal > 1
        if FMTO.FMTO_scale <= 5
            % Gain is already at minimum
            if readVal > modBoundHigh
                % Warn if power is still outside of bounds
                warning(['Femto gain is at Minimum but power is still too High\n' ...
                      ' Current Gain = ~10^%i   |  Current Power = %f'], ...
                      FMTO.FMTO_scale, readVal)
            end
            % Exit function
            return
        end
        
        % Decrease the gain
        FMTO.FMTO_scale = FMTO.FMTO_scale - 1;
        
        % Change gain
        VFN_FMTO_LUCI_setGain(FMTO);
        
        % Check power at new gain 
        read = nan(10,1);
        for i = 1:length(read)
            read(i) = FMTO.s.inputSingleScan;
        end
        readVal = mean(read(:));        
    end
    % New gain setting is done
    return
end

end
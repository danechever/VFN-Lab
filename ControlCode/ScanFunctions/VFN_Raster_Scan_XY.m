function [meas, scales, measScl, FMTO] = VFN_Raster_Scan_XY(dist, PIdevs, FMTO, delay, verbose)
% dist is a 2D matrix:
%   col1 = distX, col2 = distY
% PIdevs
% FMTO
% delay         (optional - default 0.1)
% verbose       (optional - default false)

%% Deal with input variables
%- extract desired location vectores from dist matrix
distX = dist(:,1);
distY = dist(:,2);

%- if Nread was not provided, use default
% if ~exist('Nread', 'var') || isempty(Nread)
%     Nread = 100;
% end

%- if isAutoScale was not provided, use default
% if ~exist('isAutoScale', 'var') || isempty(isAutoScale)
%     isAutoScale = true;
% end

%- if delay was not provided, use default
if ~exist('delay', 'var') || isempty(delay)
    delay = 0.1;
end

%- if verbose wasn't provided, use default
if ~exist('verbose', 'var') || isempty(verbose)
    verbose = false;
end

%% Preallocate output matrices
meas    = nan(length(distX), length(distY), FMTO.Nread);
measScl = meas;                                     % matrix for semi-processed data
scrlSz   = size(meas); sclSz(3) = 3;
scales  = nan(sclSz);                               % matrix for scales, gains, and biases

%% Perform scan
for j = 1:length(distX)
    % Check X bounds
    if (PIdevs.fibX_lower <= distX(j)) && (distX(j) <= PIdevs.fibX_upper)
        % Move fiber in x
        VFN_PIStage_move(PIdevs.fibX, distX(j));
    else
        error('fibX position out of bounds');
    end
    
    % Iterate over Y
    for i = 1:length(distY)
        % Print scan status if verbose
        if verbose
            fprintf('Scanning (%i, %i) of (%i, %i)\n', j, i, length(distX), length(distY));
        end
        
        % Check Y bounds
        if (PIdevs.fibY_lower <= distY(i)) && (distY(i) <= PIdevs.fibY_upper)
            % Move fiber in y
            VFN_PIStage_move(PIdevs.fibY, distY(i));
        else
            error('fibY position out of bounds');
        end
        
        % Let femto voltage settle at new position
        pause(delay);
        
        % Take a sample reading at current power
        read = VFN_FMTO_readN(FMTO);
        
        % Perform auto scaling if requested
        if FMTO.isAutoScale
            % Save old FMTO_scale for comparison after autoGain
            old_scl = FMTO.FMTO_scale;
            % Modify gain
            FMTO = VFN_FMTO_setAutoGain(FMTO, mean(read));
            % Check if re-read is necessary
            if old_scl ~= FMTO.FMTO_scale
                % FMTO_scale changed so the gain changed
                if verbose
                    fprintf('Femto gain was changed from %i to %i\n', old_scl, FMTO.FMTO_scale);
                end
                read = VFN_FMTO_readN(FMTO);
                if FMTO.FMTO_scale > 9
                    warning('Gain >9: %i', FMTO.FMTO_scale)
                end
            end
        end
        if find(read > 10)
            warning('Power is too high')
        end
        
        % Store data
        meas(j,i,:) = read;
        scales(j,i,1) = FMTO.FMTO_scale;
        scales(j,i,2) = FMTO.gnFact^-(FMTO.FMTO_scale-6);
        
        % Store bias and bias corrected reads
        scales(j,i,3) = VFN_FMTO_getBias(FMTO.FMTO_scale);
        measScl(j,i,:) = read - scales(j,i,3);
    end
end

end
        
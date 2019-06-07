function measScl = CourseRaster_FullScan(MDT, dist)
% NOTE: backlash is always removed from zaber motion. Even when moving forward.
% This only searches for the peak coupling and will use X,Y, and Foc
% 'dist' is a matrix whos 1st column is distX and 2nd column is distY

global s FMTO_scale

%-- Define course raster properties
distX = dist(:,1); distY = dist(:,2);
% Set backlash for piezos
backlash = 5;
% Set gainfactor for femto
gnFact  = 9.97;
% Define number of reads
Nread = s.NumberOfScans;

% Preallocate data matrices
meas=zeros(length(distX),length(distY),Nread); 
measScl = meas;         %Matrix for semi-processed data
sclSz   = size(meas); sclSz(3) = 3;
scales  = nan(sclSz);   %Matrix for scales, gains, and biases

fprintf('Starting course raster\n')
% Remove X backlash
DE2_MDTVol(MDT, distX(1)-backlash, 'x', 0);
for j=1:length(distX)
    fprintf(' Column %i of %i\n', j, length(distX))
    % Move to new X positions
    DE2_MDTVol(MDT, distX(j), 'x', 0);
    
    % Remove Y backlash
    DE2_MDTVol(MDT, distY(1)-backlash, 'y', 0);
    for i=1:length(distY)
        % Move to new Y positoin
        DE2_MDTVol(MDT, distY(i), 'y', 0);
        pause(0.1);
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
            fprintf('Femto gain was changed from %i to %i\n',old_scale, FMTO_scale)
            read    = startForeground(s);
            if FMTO_scale > 9
                warning('Gain >9: %i',FMTO_scale)
            end
        end
        
        if find(read>10)
            warning('Power is too high')
        end
        meas(j,i,:) = read;
        scales(j,i,1) = FMTO_scale;
        scales(j,i,2) = gnFact^-(FMTO_scale-6);
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
        scales(j,i,3) = locBias;
        measScl(j,i,:) = read - locBias; %Subtract bias for semi-processed matrix
    end
end
% Re-scale data to account for gain settings
measScl = mean(measScl,3).*scales(:,:,2);

end